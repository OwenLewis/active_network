! ============================================================================
! Name        : global_param_mod.f90
! Author      : Owen Lewis
! Version     :
! Copyright   : Your copyright notice
! Description : Module containing global parameters defining a fluid simulation
! ============================================================================

module global_param_mod
!PARAMETERS CRITICAL TO THE PROBLEM
    COMPLEX, PARAMETER :: J = (0,1)
    REAL, PARAMETER :: PI=DACOS(-1.D0)
	INTEGER :: Nx = 64, Ny = 64, Mb = 180
	REAL :: hx= 1./64., hy = 1./64., Xmax = 1.0, Ymax = 1.0
	REAL :: Tmax = 5., dt = 0.005, contau = 0.125
	REAL :: mu = 1, freq = 2*Pi, con_wave_num = Pi/2, Amp = 0.05, adh_offset = 0, adh_wave_num = Pi/2.
	REAL :: base_stretch_k = 10, head_multiple = 0.1
	REAL :: lame = 10, network_slip = 100, connect_stiff = 100000
	REAL :: adhesion_slip = 10, base_adh_slip = 0.000001, adhesion_stiff = 100
	REAL :: load_gamma = 0
	INTEGER :: homogeneous = 0, isotropic = 1, adhwave = 1, contractwave = 1
	INTEGER :: translate = 1
    INTEGER :: out_steps = 1
    INTEGER :: steps
    CHARACTER(len = 10) :: runname = 'testingish'

    NAMELIST /soft_params/ Tmax, dt, mu, freq, con_wave_num, Amp, base_stretch_k, head_multiple,&
                           &load_gamma, out_steps, runname, &
                           &lame, network_slip, connect_stiff, adhesion_slip, adhesion_stiff, &
                           &adh_offset, base_adh_slip, contau, homogeneous, isotropic, adhwave, &
                           &contractwave, translate, adh_wave_num, hx, hy, Xmax, Ymax, Mb, Nx, Ny



!This type defines a set of immersed boundary points, their associated tether points,
!And all associated quantities (first & second differences, etc)
!We'll put it here to make it globaly available
type im_bnd

	!These arrays define the location of immersed boundary points, associated tehter points, internal forces generated
	!by the boundary configuration, and tehter forces applied to the boundary
	REAL, DIMENSION(:,:), ALLOCATABLE :: im_bnd_loc, stretch_force, connect_force

	!These arrays store constants that define constitutive laws for forces due to stretching, and bending of the
	!boundary, as well as deviation from the associated tether points
	REAL, DIMENSION(:), ALLOCATABLE :: stretch_kp, stretch_km, gammap, gammam

	!These arrays define first difference vectors: forward, backward, and centered differences
	REAL, DIMENSION(:,:), ALLOCATABLE :: taup, taum, tau0, norm0

	!These arrays define arc length by forward, backward, and centered differences
	REAL, DIMENSION(:), ALLOCATABLE :: dsp, dsm, ds0

	!This is a generic flag associated with each point. I will use it for various purposes
	INTEGER, DIMENSION(:), ALLOCATABLE :: flag, connect_buddy

end type im_bnd


TYPE (im_bnd) :: ref_bnd

!========================================================================================================
!==============================NOW ALL THE SHIT FOR THE IMMERSED NETWORK=================================
!========================================================================================================
!========================================================================================================
!This type defines a face of an immersed structure. It includes
!identifiers for the associated network nodes, the locations
!of link centers

type network_face
!Here we  have an array of pointers to the nodes that define the face
!as well as the links between those nodes
INTEGER, DIMENSION(1:3) :: my_links, my_nodes
!This defines the center of the face
REAL, DIMENSION(1,1:2) :: face_center = 0, ref_center = 0
!This defines the area of the face
REAL :: area = 0, ref_area = 0

end type network_face


!========================================================================================================
!This type defines a node contained in an immersed structure. It includes
!identifiers for the associated faces and links that it is contained in
!as well as the location of the point
type network_node
!Here we  have an array of integers defining the indeces (within the network type)
!used to find the faces and links associated with this node
INTEGER, DIMENSION(:), ALLOCATABLE :: my_links, my_faces
!This defines the location of the node
REAL, DIMENSION(1,1:2) :: location, elastic_force, adh_force = 0, connect_force = 0
LOGICAL :: boundary = .false., linked = .false.
REAL :: area = 0, ref_area = 0, connect_k = 0, node_drag
INTEGER :: link_count = 0, face_count = 0, connect_buddy
end type network_node

!========================================================================================================
!This type defines a link contained in an immersed structure. It includes
!identifiers for the faces it demarcates and the nodes @ either end,
!as well as some parameters related to the constitutive law of the network

type network_link
!Here we  have arrays of integers defining the indeces (within the network type)
!used to find the associated faces and nodes of this link
INTEGER, DIMENSION(1:2) :: my_faces, my_nodes
!NOTE RIGHT HERE: BE VERY AWARE THAT A LINK WITH ONE OF ITS FACES AS ZERO MEANS IT IS AN EXTERIOR LINK
!THIS IS TO SAY THAT THE ZERO'TH FACE IS THE FACE 'OUTSIDE' OF THE ENTIRE NETWORK
REAL :: link_length = 0, ref_length = 0, active_tens = 0
REAL, DIMENSION(1,1:2) :: rayup, raydown
LOGICAL :: boundary = .false.
REAL :: stiff = 0
!This defines the center of the link
REAL, DIMENSION(1,1:2) :: link_center = 0, ref_center = 0
end type network_link

!========================================================================================================
!This type defines a set of immersed network points, their associated faces, etc
!We'll put it here to make it globaly available
type im_network
	!An array of reals that defines the location of nodes in the network
	TYPE(network_node), DIMENSION(:), ALLOCATABLE :: nodes
	!An array of associated faces of the network
	type (network_face), DIMENSION(:), ALLOCATABLE :: faces
	!An array of links between nodes (aka edges of faces)
	type (network_link), DIMENSION(:), ALLOCATABLE :: links
	!Integers that define the length of these various lists
	INTEGER :: M_nodes, M_faces, M_links
	!A real jus to store the total elastic strain energy associated with the deformation of the network
	REAL :: strain_energy

end type im_network


!========================================================================================================
!========================================================================================================
!========================================================================================================
!========================================================================================================

!This type defines a set of adhesion points, their associated faces, etc
!We'll put it here to make it globaly available
type adh_complex
	!An array of reals that defines the location of nodes in the network
	REAL, DIMENSION(:), ALLOCATABLE :: adh_length, adh_stiff, adh_slip
	REAL, DIMENSION(:,:), ALLOCATABLE :: location, ref_location, adh_dir
	!An array of integers that determines which network point each adhesion is linked to
	INTEGER, DIMENSION(:), ALLOCATABLE :: network_buddy
	!Integers that define the length of these various lists
	INTEGER :: M_adh
	!And a real just to store the total strain energy associated with adhesion links
	REAL :: strain_energy

end type adh_complex


end module global_param_mod
