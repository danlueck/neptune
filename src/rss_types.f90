 module rss_types
    use slam_types,         only: sp,dp,i8b
    use slam_time,          only: time_t
    use slam_orbit_types,   only: state_t,covariance_t
    use maneuvers,          only: maneuver_t
    use neptuneClass,       only: Neptune_class

    implicit none

    integer,parameter                       :: MAX_CHEBYSHEV_DEGREE = 36        ! Maximum chebyshev degree

    type t_ephemeris
        integer(i8b)                        :: ephemeris_id                     ! ID of the ephemeris in the database

        type(time_t)                        :: epoch                            ! Epoch the ephemeris is valid for              [mjd]
        !* Keplerian elements
        real(dp)                            :: semi_major_axis  = 0.0d0         ! The semi-major axis                           [km]
        real(dp)                            :: eccentricity     = 0.0d0         ! The eccentricity                              [-]
        real(dp)                            :: inclination      = 0.0d0         ! The inclination                               [rad]
        real(dp)                            :: raan             = 0.0d0         ! The right ascension of teh ascending node     [rad]
        real(dp)                            :: aop              = 0.0d0         ! The argument of perigee                       [rad]
        real(dp)                            :: true_anomaly     = 0.0d0         ! The mean anomaly                              [rad]
        real(dp)                            :: mean_anomaly     = 0.0d0         ! The mean anomaly                              [rad]
        real(dp)                            :: ballistic_coefficient = 0.0d0    ! The ballistic coefficient                     [m^2/kg]
        !* Cartesian coordinates   
        real(dp),dimension(6)               :: state_vector     = 0.0d0         ! The x,y,z position vector                     [km]
                                                                                !  and the x_dot,y_dot,z_dot vector             [km/s]
        type(covariance_t)                  :: covariance_matrix                ! The error covariances
        logical                             :: propagated_flag  = .false.       ! Flag indicating that the ephemeris have been extrapolated
                                                                                !  and are not based on measurement values at epoch
!        type(t_propagator)                  :: propagator                       ! Propagator used for extrapolation
        integer,dimension(:,:),allocatable      :: polynomial_degree            ! Degree of the chebyshev polynomials
        real(dp),dimension(:,:,:,:),allocatable :: polynomial_coefficients      ! Used for interpolation (e.g. Chebychev)       [-]
        type(time_t),dimension(:,:),allocatable :: interpolation_start_epoch    ! The starting point of interpolation for each granule
        type(time_t),dimension(:,:),allocatable :: interpolation_end_epoch      ! The end point of interpolation for each granule
        logical                                 :: is_manoeuvering  = .false.   ! Indicates that the satellite is manoeuvering between 
                                                                                !  interpolation start and end epoch.
    end type t_ephemeris

    ! Attributes of a manoeuvre
    type t_manoeuvre
        type(time_t)                        :: start_epoch                      ! Date of maneuver start
        type(time_t)                        :: end_epoch                        ! Date of maneuver end (will be determined based on duration or thrust ratio)
        logical                             :: has_thrust_ratio                 ! Indicates that a thrust ratio rather than a duration is used
        real(dp)                            :: duration            = 0.0d0      ! Duration of manoeuvre                         [s]
        integer                             :: no_thrust_revs      = 0          ! Revolutions with no-thrust                    [-]
        integer                             :: thrust_revs         = 0          ! Revolutions with thrust                       [-]
        real(dp),dimension(3)               :: acc                              ! Acceleration vector in UVW frame              [m/s**2]
        real(dp)                            :: rel_uncertainty                  ! Relative uncertainty per ignition             [-]
    end type t_manoeuvre

        type :: t_pass
        integer                             :: object_id      = 0               ! ID of object
        type(time_t)                        :: epoch                            ! Epoch the ephemeris is valid for  [mjd]
        real(dp),dimension(6)               :: rae_vector     = 0.0d0           ! The RAE position vector           [km],[rad],[rad],[km/s],[rad/s],[rad/s]
    end type t_pass

    ! Attributes for object
    type :: t_object
        integer                             :: object_id                        ! ID that identifies the link to the catlog
        integer                             :: norad_id                         ! Object identified by NORAD id
        character(len=255)                  :: name                             ! Defines which population this object belongs to
        character(len=255)                  :: population_identifier            ! Defines which population this object belongs to
        integer                             :: numEphemerides    = 0            ! Number of ephemerides available
        integer                             :: numManoeuvres     = 0            ! Number of manoeuvres available
        integer                             :: numPasses         = 0            ! Number of passes available

        logical                             :: validEphemerides  = .false.      ! Flag that indicates whether the OD worked
                                                                                !  worked successfully.
        integer                             :: validEphemerisIndex = 0          ! Flag that indicates whether the OD worked
                                                                                !  worked successfully.
        logical                             :: reentered_flag    = .false.      ! Indicates that an object reentered into the Earth's atmosphere
        real(dp)                            :: mass              = 8110.0d0     ! Mass of object                                [kg]
        real(dp)                            :: diameter          = 9.732d0      ! Diameter of object                            [m]
        real(dp)                            :: rcs               = 74.39d0      ! Radar cross section                           [m^2]
        real(dp)                            :: area_to_mass      = 9.1726d-3    ! Area to mass ratio                            [m^2/kg]
        real(dp)                            :: c_D               = 2.2d0        ! Drag coefficient                              [-]
        real(dp)                            :: c_R               = 1.3d0        ! Reflectivity coefficient                      [-]
        real(dp)                            :: factor            = 1.0d0        ! Factor ID from MASTER population.
        type(t_ephemeris),dimension(:),allocatable  :: ephemerides              ! The ephemerides of the object
        type(t_manoeuvre),dimension(:),allocatable  :: manoeuvres               ! Manoeuvres of the object
        integer                             :: last_mnv_idx      = 1            ! Index of last considered manoeuvre
        real(dp)                            :: remaining_duration= 0.0d0        ! For a thrust-ratio manoeuvre that may span multiple 
                                                                                !  propagation steps the remaining duration is stored.
        logical                             :: thrust_phase      = .false.      ! Indicates in which phase pf the thrust-ratio manoeuvre we left off.
        type(maneuver_t),dimension(:),allocatable   :: neptune_manoeuvres       ! Manoeuvre data srtucture as used in NEPTUNE internally

        ! CAMP specific
        type(t_pass),dimension(:),allocatable       :: passes                   ! Passes expressed as RAE vector for wrt a given sensor
                                                                                !  [km],[rad],[rad],[km/s],[rad/s],[rad/s]
        character(len=20)                   :: timestamp                        ! Time stamp later used to identify the passes file
        integer                             :: missed_interpolations = 0        ! Counts the number of missed interpolations attempts

        ! MWG specific
        logical                             :: wrongHemisphere   = .false.      ! Indicates that the satellite is no where near the FOV
        real(dp)                            :: hemisphereTimeout = 0.0d0        ! For optimization purposes we keep track how long we are skipping the propagation
        integer                             :: exclusive_sensor  = -1           ! Sensor ID, which identifies the senors which has the exclusive detection rights on this object (e.g. GPS sensor enabled, no need for additional RADAR or telescope measurment)
    end type t_object

    ! Attributes for radar detection
    type t_detection
        integer                             :: object_id       = 0              ! ID of object
        integer                             :: detection_id    = 0              ! ID of detection
        integer                             :: tracklet_id     = 0              ! ID of corresponding tracklet
        integer                             :: tracking_switch = 0              ! 0/1 = notrack/track
        integer                             :: tracklet_score  = 0              ! Score of detection in tracklet
        type(time_t)                        :: time                             ! Time in Julian Days
        logical                             :: detected        = .false.        ! Defines whether an objects has been detected
        real(dp),dimension(6)               :: state_vector    = 0.0d0          ! x,y,z,dx,dy,dz in inertial    [km,km/s]
        real(dp),dimension(6)               :: ecef_state      = 0.0d0          ! x,y,z,dx,dy,dz in earth-fixed [km,km/s]
        real(dp)                            :: azimuth         = 0.0d0          ! Azimuth                       [rad]
        real(dp)                            :: elevation       = 0.0d0          ! Elevation                     [rad]
        real(dp)                            :: range           = 0.0d0          ! Range                         [km]
        real(dp)                            :: range_rate      = 0.0d0          ! Range rate                    [km/s]
        real(dp)                            :: azimuth_rate    = 0.0d0          ! Azimuth rate                  [rad/s]
        real(dp)                            :: elevation_rate  = 0.0d0          ! Elevation rate                [rad/s]
        integer                             :: num_integrations= 0              ! # of integratinos             [-]
        real(dp)                            :: snr             = 0.0d0          ! Signal to noise ration        [dB]
        real(dp)                            :: rcs             = 0.0d0          ! Radar cross section           [dBsm]
        real(dp)                            :: probability     = 0.0d0          ! Probability of detection      [-]
        real(dp)                            :: ballistic_coefficient = 0.0d0    ! ballistic coefficient = c_D * A/m [m^2/kg]
        logical                             :: is_manoeuvering  = .false.       ! Indicates that the satellite is manoeuvering
    end type t_detection

    ! Attributes for radar tracklet
    type t_tracklet
        !sequence
        integer                                :: tracklet_id = 0               ! ID for tracklet, see detection
        integer                                :: num_detections  = 0           ! Number of detections forming the tracklet
        type (t_detection),dimension(:),allocatable :: tracklet_detections      ! Detections in tracklet
    end type t_tracklet

    ! Attributes for radar station
    type :: t_station

        character(len=255)              :: name                                 ! Name of the station
        integer                         :: antenna_type     = 0                 ! Type of Antenna (0 = Parabol, 
                                                                                !                   1 = Phased-Array, 
                                                                                !                   2 = Telescope, 
                                                                                !                   3 = Laser Ranging, 
                                                                                !                   4 = GPS antenna)

        integer                         :: operation_mode   = 0                 ! Operation mode: 0 - Basic Mode (IRAS implementation)
                                                                                !                 1 - Mechanical Tracking Mode (*.mtr) - FHR OVER
                                                                                !                 2 - Electronic Tracking Mode (*.etr) - FHR OVER
                                                                                !                 3 - Scanning Mode            (*.scn) - FHR OVER
                                                                                !                 4 - Optical Tracking         (*.otr) - IRAS OPM
                                                                                !                 5 - GPS sensor               (*.gps) - IRAS GPS
        character(len=255)              :: operation_mode_name
        logical                         :: noise_flag      = .false.            ! Flag indicating the use of noisy observations
        real(dp),dimension(6)           :: measurement_noise                    ! Diagonals of the measruement noise matrix in [km],[km/s] and [rad/s]
        logical                         :: snr_noise_flag  = .false.            ! Flag indicates the use of the range, range rate and angle sigma
                                                                                !  values instead of the measurement noise values.
        real(dp)                        :: altitude     = 0.0d0                 ! Altitude of station          [m]
        real(dp)                        :: longitude    = 0.0d0                 ! Longitude of station         [rad]
        real(dp)                        :: latitude     = 0.0d0                 ! Latitude of station          [rad]
        real(dp)                        :: max_range    = 0.0d0                 ! Maximum detection range      [km]

        ! Available after initialization through site2ecef(site) subroutine
        real(dp),dimension(3)           :: r_ecef       = 0.0d0                 ! Position vector of the site  [km]
                                                                                !  in the earth centered earth
                                                                                !  fixed (ECEF) frame.
        !* Cartesian coordinates
        real(dp),dimension(6)           :: state_vector     = 0.0d0             ! The x,y,z position vector                     [km]
                                                                                !  and the x_dot,y_dot,z_dot vector             [km/s]
        real(dp)                        :: rightAscension0                      ! Sensor's right ascension at initilization [rad]
        real(dp)                        :: rightAscension                       ! Sensor's right ascension at epoch  [rad]
        
        ! Helper variables
        real(dp)                    :: SIN_LAT          = 0.0d0                 ! SIN value of the latitude angle  [rad]
        real(dp)                    :: COS_LAT          = 0.0d0                 ! COS value of the latitude angle  [rad]
        real(dp)                    :: SIN_LONG         = 0.0d0                 ! SIN value of the longitude angle [rad]
        real(dp)                    :: COS_LONG         = 0.0d0                 ! COS value of the longitude angle [rad]

        real(dp)                    :: C_sensor                                 ! location specific constant which makes later updates in ECI much faster

        type(t_detection),dimension(:),allocatable  :: detections               ! Array of detections
        integer                                     :: num_detections = 0       ! Counts the # of detcetions in the detections array
        integer                                     :: max_detections = 1000000 ! Defines the size of the detections array
        type(t_tracklet),dimension(:),allocatable   :: tracklets                ! Array of tracklets
        integer                                     :: num_tracklets = 0        ! describes how many tracklets are in the array at a given time
        integer                                     :: all_tracklets = 0        ! overall tracklet statistics
        integer                                     :: all_reliable_tracklets = 0


    end type t_station

        ! Attributes of a conjunction
    type :: t_conjunction
        character(len=255)          :: analysis_details                         ! Details of the analysis
        type(time_t)                :: epoch                                    ! Time of close approach (tca)
        type(time_t)                :: report_epoch                             ! Time of report generation
        real(dp)                    :: tca_std                 = 0.d0           ! 1-sigma std in time [s]
        type(t_object),pointer      :: object_1                                 ! Object 1 of conjunction pair
        type(t_object),pointer      :: object_2                                 ! Object 2 of conjunction pair
        !type(t_ephemeris)           :: ephemeris_1                              ! Ephemeris 1
        !type(t_ephemeris)           :: ephemeris_2                              ! Ephemeris 2
        integer(i8b)                :: ephemeris_id_1           = -1            ! Ephemeris ID of target object as used in the database
        integer(i8b)                :: ephemeris_id_2           = -1            ! Ephemeris ID of risk object as used in the database
        integer                     :: object_id_1              = -1            ! Object ID of target as used in the database
        integer                     :: norad_id_1               = -1            ! Norad ID of target
        character(len=255)          :: object_type_1                            ! Type of the target object (e.g. PAYLOAD, DEBRIS, ...)
        integer                     :: object_id_2              = -1            ! Object ID of risk object as used in the database
        integer                     :: norad_id_2               = -1            ! Norad ID of risk object
        character(len=255)          :: object_type_2                            ! Type of the risk object (e.g. PAYLOAD, DEBRIS, ...)
        real(dp)                    :: diameter_1               = 0.0d0         ! Diameter of target object [m]
        real(dp)                    :: diameter_2               = 0.0d0         ! Diameter of risk object [m]
        real(dp)                    :: cross_section_1          = 0.0d0         ! Cross section of target object [m**2]
        real(dp)                    :: cross_section_2          = 0.0d0         ! Cross section of risk object [m**2]
        real(dp),dimension(6)       :: state_vector_1           = 0.0d0         ! Cartesian state vector of target
        real(dp),dimension(6)       :: state_vector_2           = 0.0d0         ! Cartesian state vector of risk object
        real(dp),dimension(6)       :: keplerian_elements_1     = 0.0d0         ! Keplerian elements of target
        real(dp),dimension(6)       :: keplerian_elements_2     = 0.0d0         ! Keplerian elements of risk object
        real(dp),dimension(3)       :: covariance_elem_1        = 0.0d0         ! Covariance elements of target
        real(dp),dimension(3)       :: covariance_elem_2        = 0.0d0         ! Covariance elements of risk object
        real(dp)                    :: close_approach_k_c_SQ    = -1.d0
        real(dp),dimension(6)       :: delta_vector_eci         =  0.d0         ! Delta vector of objects in ECI    [km] and [km/s]
        real(dp),dimension(6)       :: delta_vector_uvw         =  0.d0         ! Delta vector of objects in UVW    [km] and [km/s]
        real(dp),dimension(6)       :: delta_vector_uvw_std     =  0.d0         ! Std of delta vector of objects in UVW    [km] and [km/s]
        real(dp)                    :: p_col_max                = 0.d0          ! Maximum collision probability     [-]
        real(dp)                    :: p_col_mc                 = 0.d0          ! Monte Carlo collision probability [-]
        real(dp)                    :: p_col_alfriend           = 0.d0          ! Alfriend and Akella collision probability     [-]
        real(dp)                    :: p_col_chan               = 0.d0          ! Chan collision probability     [-]
        real(dp)                    :: p_col_foster             = 0.d0          ! Foster collision probability     [-]
        real(dp)                    :: p_col_patera             = 0.d0          ! Patera collision probability     [-]
        real(dp)                    :: p_col_alfano             = 0.d0          ! Alfano collision probability     [-]
    end type t_conjunction

     ! Attributes of a conjunctions report
    type :: t_conjunctions
        integer                     :: report_id                = -1            ! Identifies the report in the database
        type(time_t)                :: report_epoch                             ! Indicates the report epoch
        character(:),allocatable    :: target_population_id                     ! Identifies the population of target objects
        character(:),allocatable    :: risk_population_id                       ! Identifies the population of risk objects
        integer                     :: num_target_objects       = 0             ! Number of target objects screened
        integer                     :: num_risk_objects         = 0             ! Number of risk objects screened
        integer                     :: num_filtered_risk_objects= 0             ! Number of filtered risk objects
        integer                     :: num_conjunctions         = 0             ! Number of conjunctinos in the list
        type(t_conjunction),dimension(:),allocatable    :: list                 ! List of conjunctions

    end type t_conjunctions

end module rss_types
