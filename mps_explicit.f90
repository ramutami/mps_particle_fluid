!------------change allowed--------------!
    !module parameters_for_simulationを適宜変更すること。
    !module inital_condition下に初期条件について定めるsubroutineを作る必要がある。
    !*all global parameters and variables must be defind in module:parameters_and_variables_for_simulation
!------------change allowed--------------!

!------------Caution!--------------!
    !i,j,j…などはsubroutine内部でのみ使うようにすること。
    !３次元で計算する場合は重力加速度を方向に設定すること。
!------------Caution!--------------!


module mpi_module
    !mpiに関するmodule
    use mpi_f08
    implicit none

    integer :: myrank, nprocs
    integer :: i_start, i_finish, n_local
    integer, allocatable :: counts(:), displs(:)

    contains

    subroutine setup_mpi_decomposition(n)
        implicit none
        integer, intent(in) :: n
        integer :: ierr
        integer :: r, base, rem, offset

        call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

        allocate(counts(nprocs))
        allocate(displs(nprocs))

        base = n / nprocs
        rem  = mod(n, nprocs)

        offset = 0
        do r = 0, nprocs-1
            counts(r+1) = base
            if (r < rem) counts(r+1) = counts(r+1) + 1

            displs(r+1) = offset
            offset = offset + counts(r+1)
        end do

        n_local = counts(myrank+1)
        i_start = displs(myrank+1) + 1
        i_finish   = i_start + n_local - 1

    end subroutine setup_mpi_decomposition

    subroutine allgather_real_vector(a)
        implicit none
        real(8), intent(inout) :: a(:)
        real(8), allocatable :: sendbuf(:)
        integer :: ierr

        allocate(sendbuf(n_local))
        sendbuf(:) = a(i_start:i_finish)

        call MPI_Allgatherv( &
            sendbuf, n_local, MPI_DOUBLE_PRECISION, &
            a, counts, displs, MPI_DOUBLE_PRECISION, &
            MPI_COMM_WORLD, ierr)

        deallocate(sendbuf)
    end subroutine allgather_real_vector


    subroutine allgather_integer_vector(a)
        implicit none
        integer, intent(inout) :: a(:)
        integer, allocatable :: sendbuf(:)
        integer :: ierr

        allocate(sendbuf(n_local))
        sendbuf(:) = a(i_start:i_finish)

        call MPI_Allgatherv( &
            sendbuf, n_local, MPI_INTEGER, &
            a, counts, displs, MPI_INTEGER, &
            MPI_COMM_WORLD, ierr)

        deallocate(sendbuf)
    end subroutine allgather_integer_vector

    subroutine allgather_real_vector3(a)
        implicit none
        real(8), intent(inout) :: a(:,:)
        real(8), allocatable :: sendbuf(:)
        integer :: ierr, k

        allocate(sendbuf(n_local))

        do k = 1, 3
            sendbuf(:) = a(i_start:i_finish,k)

            call MPI_Allgatherv( &
                sendbuf, n_local, MPI_DOUBLE_PRECISION, &
                a(1,k), counts, displs, MPI_DOUBLE_PRECISION, &
                MPI_COMM_WORLD, ierr)
        end do

        deallocate(sendbuf)
    end subroutine allgather_real_vector3

end module mpi_module

module parameters_and_variables_for_simulation
    implicit none


    !gloval_parameters
        !dimention: シミュレーションの次元を選択する。2->二次元、３->三次元
        !pargicle_distance: 初期の粒子間距離[m]
        !time_interval: 計算ステップの時間幅[s]、一ステップでは2*time_intervalだけ時間が進む。
        !output_interval: 何ステップごとに粒子の位置を出力するか
        !finish_time: 終了時刻[s]
        !viscosity:　粘性
        !collision_distance_ratio: 衝突判定のさいの粒子距離の関値に関するパラメタcollision_distance/particle_distance
        !restitution_coefficient: 反発係数
        !Re_for_..._parameter: 重み関数で使う臨界距離Re_for_laplacian,Re_for_numberdensity,etcに対し、Re_for..=Re_for_..._parameter*particle_distance
        !threshold_ratio_of_number_density: number_density< (threshh...) * n0_for_number_densityの時に表面粒子と判定する。
        !relaxation_coefficient_for_pressure：ポアソン方程式のソース項の修正係数。
        !compressibility：ポアソン方程式の行列の修正項の大きさ。正定値性を強めてくれる。
        !fluid_density：計算で使う流体密度rho0の値。
        !cg_max_iterr：共役勾配法で解く際のmaxiteration.
        !cg_relative_tolerance：共役勾配法のbreak条件
        !umax：想定される最大流速

        !------------change allowed--------------!
        integer :: dimention = 3
        real(8),parameter :: particle_distance=0.025
        real(8),parameter :: time_interval=0.001
        integer,parameter :: max_timestep = 20000
        integer, parameter :: output_interval=20
        real(8),parameter :: finish_time = 2.0
        real(8),parameter :: viscosity = 1.0e-6
        real(8), parameter :: collision_distance_ratio = 0.5
        real(8),parameter :: restitution_coefficient = 0.2
        real(8),parameter :: fluid_density = 1000.0_8
        real(8),parameter :: Re_for_laplacian_parameter = 3.1
        real(8),parameter :: Re_for_number_density_parameter = 2.1
        real(8),parameter :: Re_for_gradient_parameter = 2.1
        real(8),parameter :: threshold_ratio_of_number_density = 0.97
        real(8),parameter :: relaxation_coefficient_for_pressure = 0.2
        real(8),parameter :: compressibility = 0.45e-9
        integer,parameter :: cg_max_iterr = 500
        real(8),parameter :: cg_relative_tolerance = 1.0e-6
        real(8),parameter :: umax = 4.9
        !------------change allowed--------------!


        integer,parameter :: fluid = 1
        integer,parameter :: wall = 2
        integer,parameter :: dummywall = 3
        integer,parameter :: surface_particle = 1
        integer,parameter :: inner_particle = 0
        integer,parameter :: dummy = -1
        real(8),parameter :: g_x = 0.0
        real(8),parameter :: g_y = -9.80665
        real(8),parameter :: g_z = 0.0
        real(8),parameter :: Re_for_laplacian = Re_for_laplacian_parameter*particle_distance
        real(8),parameter :: Re_for_gradient = Re_for_gradient_parameter*particle_distance
        real(8),parameter :: Re_for_number_density = Re_for_number_density_parameter * particle_distance
        real(8),parameter :: collision_distance = collision_distance_ratio*particle_distance
        real(8),parameter :: sound_speed_for_calculation = 5.0*umax
    !

    !global_variables
        !global変数：すべてのmoduleのsubroutineで参照可能な変数
        !allocatはsubroutineで行う。

        !particle_position: 粒子iのx,y,z座標をparticle_position(i,x=1/y=2/z=3)に収納する。
        !particle_velocity: 粒子iのx,y,z方向の速度ををparticle_velocity(i,x=1/y=2/z=3)に収納する。
        !particle_acceleration: 粒子iに関する加速度をacceleration(i,x=1/y=2/z=3)に収納する。
        !number_density: 重みつき粒子数密度Σw(ri-rj)
        !boundary_condition: 自由表面境界の判別
        !number_of_particles：粒子総数。粒子の初期配置を決定した段階で、particle_position等はこの数にallocateする。
        !particle_type：粒子の属性。液体粒子か、壁粒子か、ダミー壁粒子か
        !particle_prssure：各粒子位置での圧力
        !original_layer：各粒子が最初、水柱の中でどの高さに属しているか。
        !source_term

        real(8),allocatable :: particle_position(:,:)
        real(8),allocatable :: particle_velocity(:,:)
        real(8), allocatable :: velocity_after_collision(:,:)
        real(8),allocatable :: particle_acceleration(:,:)
        real(8),allocatable :: number_density(:)
        integer,allocatable :: boundary_condition(:)
        integer,allocatable :: particle_type(:)
        real(8),allocatable :: particle_pressure(:)
        real(8),allocatable :: Original_layer(:)
        integer :: number_of_particles
        real(8) :: n0_for_laplacian
        real(8) :: n0_for_number_density
        real(8) :: n0_for_gradient
        real(8) :: lambda_0
        real(8),allocatable :: source_term(:)
        real(8),allocatable :: minimum_pressure(:)

    !

end module parameters_and_variables_for_simulation

module initial_particle_position_velocity_particle_type
    use omp_lib
    use parameters_and_variables_for_simulation
    implicit none

    !以下のsubroutineで初期の粒子の配置、壁の配置、ダミー壁の配置、粒子速度を定める。
    !壁、ダミー壁はすべて粒子として扱う。
    !粒子についての位置particle_positions(dim,n)を定めた上で、
    !その粒子のtype:particle_type(n)が液体粒子なのか壁なのかダミー壁なのかを(後から)定めることで、
    !液体粒子、壁、ダミー壁の位置を定めることができる。

    contains
    subroutine water_tank_and_water_column_2d(x_watertank,y_watertank,x_watercolumn,y_watercolumn,wallthickness,dummywallthickness)

        !二次元水柱崩壊を計算する場合の初期条件
        !水槽の大きさは,壁の内側では買った大きさがx_watertank*y_watertank、水柱の大きさはx_watercolumn*y_watercolumn [m]

        use parameters_and_variables_for_simulation
        implicit none

        real(8),intent(in)::x_watertank,y_watertank,x_watercolumn,y_watercolumn
        integer,intent(in)::wallthickness,dummywallthickness
        !水槽や水柱のx,y方向の粒子数、例えばnx_watertank = x_watertank/particle_distance
        integer :: nx_watertank,ny_watertank,nx_watercolumn,ny_watercolumn
        real(8) :: epsilon = 0.01
        !計算ようの変数
        integer :: i,ix,iy
        integer :: temp_int_1,temp_int_2,temp_int_3

        !nx,nyの計算。epsilon*particl_distanceをつけるのは、real(8)の不安定さ故。
        nx_watertank = floor((x_watertank+particle_distance*epsilon)/particle_distance)
        ny_watertank = floor((y_watertank+particle_distance*epsilon)/particle_distance)
        nx_watercolumn = floor((x_watercolumn+particle_distance*epsilon)/particle_distance)
        ny_watercolumn = floor((y_watercolumn+particle_distance*epsilon)/particle_distance)
        write(*,*) "water tank size :", nx_watertank,ny_watertank
        write(*,*) "water column size : ",nx_watercolumn,ny_watercolumn

        !num_of_particle
        !必要な粒子数のメモリの数：
        !壁については　(nxtank+2*(wallthickness+dummywallthickness))*(nytank+wallthicnkness+dummywallthickness)-nx_tank*ny_tank
        !水柱については　nx_watercolumn*ny_watercolumn
        number_of_particles = (nx_watertank+2*(wallthickness+dummywallthickness))*(ny_watertank+wallthickness+dummywallthickness)-&
        & nx_watertank*ny_watertank
        number_of_particles = number_of_particles + nx_watercolumn*ny_watercolumn
        write(*,*) 'number of particles : ', number_of_particles

        !メモリの確保
        !二次元三次元によらず、三次元の座標を確保する。（計算の内容は変わらないため。）
        allocate(particle_position(number_of_particles,3))
        allocate(particle_velocity(number_of_particles,3))
        allocate(particle_acceleration(number_of_particles,3))
        allocate(velocity_after_collision(number_of_particles,3))
        allocate(number_density(number_of_particles))
        allocate(boundary_condition(number_of_particles))
        allocate(particle_type(number_of_particles))
        allocate(particle_pressure(number_of_particles))
        allocate(Original_layer(number_of_particles))
        allocate(source_term(number_of_particles))
        allocate(minimum_pressure(number_of_particles))

        velocity_after_collision(:,:) = 0.0
        number_density(:) = 0.0_8
        boundary_condition(:) = 0
        particle_pressure(:) = 0.0
        minimum_pressure(:) = 0.0

        !初期粒子位置の設定        
        !y座標を固定し、xについてのloopを回しながら初期条件を設定してゆく。

        !以下はdummywall＋wallの厚さが4の場合

        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |○____________________________________________________|           | iY=+1
        !|                                                                              | iY=-0
        !|                                                                              | iY=-1
        !|                                                                              | iY=-2 
        !|______________________________________________________________________________| iY=-3
        ! -3,-2,-1, 0,+1 =iX

        !○の位置が(ix,iy) = 1,1の座標となる。
        !粒子番号iについては、上の図で左下を1として、
        !xを走査しながらyを上げていってiを順次定めていく。（water_tankに対し）
        !また、water_tankの最後の粒子をi=i_last_tankとした時、
        !水柱の左下に来る粒子をi=i_last_tank+1として、同じようにxを走査しながらyを上げていってiを順次定めていく。

        !計算量削減のため、毎回計算する変数はループの外側でtemp_intで計算してしまう。

        !床dummy
        temp_int_1 = 1-wallthickness-dummywallthickness
        temp_int_2 = nx_watertank+2*(wallthickness+dummywallthickness)
        temp_int_3 = 1-wallthickness-dummywallthickness
        !$omp parallel do private(i,iX)
        do iY = 1-wallthickness-dummywallthickness,1-wallthickness-1
            do iX = 1-wallthickness-dummywallthickness,nx_watertank+wallthickness+dummywallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do 
        end do
        !end $omp parallel do

        !床wall/dummy
        !$omp parallel do private(i,iX)
        do iY = 1-wallthickness,1-1
            do iX = 1-wallthickness-dummywallthickness,1-wallthickness-1    
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = 1-wallthickness,nx_watertank+wallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+wallthickness+1,nx_watertank+wallthickness+dummywallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do

        !壁
        temp_int_1 = (wallthickness+dummywallthickness)*(nx_watertank+2*(wallthickness+dummywallthickness))
        temp_int_2 = 2*(wallthickness+dummywallthickness)
        temp_int_3 = 1-wallthickness-dummywallthickness
        !$omp parallel do private(i,iX)
        do iY = 1,ny_watertank-wallthickness
            do iX = 1-wallthickness-dummywallthickness,1-wallthickness-1
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = 1-wallthickness,1-1
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+1,nx_watertank+wallthickness
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+wallthickness+1,nx_watertank+wallthickness+dummywallthickness
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do
        !$omp parallel do private(i,iX)
        do iY = ny_watertank-wallthickness+1,ny_watertank
            do iX = 1-wallthickness-dummywallthickness,1-1
            i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = wall
            Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+1,nx_watertank+wallthickness+dummywallthickness
            i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = wall
            Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do
        
        temp_int_1 = (nx_watertank+2*(wallthickness+dummywallthickness))*(ny_watertank+wallthickness+dummywallthickness)-&
        & nx_watertank*ny_watertank
        do iY = 1,ny_watercolumn   
            do iX = 1,nx_watercolumn
            i = temp_int_1 + (iY-1)*nx_watercolumn+iX
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = fluid
            Original_layer(i) = iY*particle_distance
            end do
        end do

    end subroutine


    subroutine water_tank_and_water_column_3d(x_watertank,y_watertank,z_watertank,x_watercolumn_min,y_watercolumn_min,z_watercolumn_min,x_watercolumn_max,y_watercolumn_max,z_watercolumn_max,wallthickness,dummywallthickness)
        use parameters_and_variables_for_simulation
        implicit none

        !水槽の内側は(0,0,0)→(x_watertank,y_watertank,z_watertank)
        !水柱は(x_watercolumn_min,y_watercolum_min,z_watercolumn_min)→(x_watercolumn_max,y_watercolum_max,z_watercolumn_max)
        !水柱が水槽（＋その上空）に収まるように指定する必要がある。
        real(8),intent(in)::x_watertank,y_watertank,z_watertank,x_watercolumn_min,y_watercolumn_min,z_watercolumn_min,x_watercolumn_max,y_watercolumn_max,z_watercolumn_max
        integer,intent(in)::wallthickness,dummywallthickness

        !水槽や水柱のx方向の粒子数、例えばnx_watertank = x_watertank/particle_distance
        integer :: nx_watertank,ny_watertank,nz_watertank,nx_watercolumn,ny_watercolumn,nz_watercolumn
        real(8) :: epsilon = 0.01

        !計算用の変数
        integer :: i
        real(8) :: x,y,z
        real(8) :: max_z

        !nx,nyの計算。epsilon*particl_distanceをつけるのは、real(8)の不安定さ故。
        nx_watertank = floor((x_watertank+particle_distance*epsilon)/particle_distance)
        ny_watertank = floor((y_watertank+particle_distance*epsilon)/particle_distance)
        nz_watertank = floor((z_watertank+particle_distance*epsilon)/particle_distance)
        nx_watercolumn = floor((x_watercolumn_max-x_watercolumn_min+particle_distance*epsilon)/particle_distance)
        ny_watercolumn = floor((y_watercolumn_max-y_watercolumn_min+particle_distance*epsilon)/particle_distance)
        nz_watercolumn = floor((z_watercolumn_max-z_watercolumn_min+particle_distance*epsilon)/particle_distance)
        write(*,*) "water tank size :", nx_watertank,ny_watertank,nz_watertank
        write(*,*) "water column size : ",nx_watercolumn,ny_watercolumn,nz_watercolumn

        !num_of_particle:必要な粒子数のメモリの数
        number_of_particles = 0
        max_z = max(z_watertank,z_watercolumn_max)
        z = 0.0-real((wallthickness+dummywallthickness),8)*particle_distance
        do while (z<=max_z+particle_distance*epsilon)
            y = 0.0-(wallthickness+dummywallthickness)*particle_distance

            do while (y<=(y_watertank+real((wallthickness+dummywallthickness),8)*particle_distance)+particle_distance*epsilon)

                x = 0.0-(wallthickness+dummywallthickness)*particle_distance
                do while (x<=x_watertank+(real((wallthickness+dummywallthickness),8)*particle_distance)+particle_distance*epsilon)

                    if (z<=z_watertank .and. (x<0 .or. y<0 .or. z<0 .or. x>x_watertank .or. y>y_watertank) ) then
                        !床・壁の粒子数
                        number_of_particles = number_of_particles + 1
                    else if (x>=x_watercolumn_min .and. y>=y_watercolumn_min .and. z>=z_watercolumn_min .and. x<= x_watercolumn_max .and. y<=y_watercolumn_max .and. z<= z_watercolumn_max) then
                        !流体粒子
                        number_of_particles = number_of_particles + 1
                    end if

                    x = x + particle_distance
                end do
                y = y + particle_distance
            end do
            z = z + particle_distance
        end do

        write(*,*) 'number of particles : ', number_of_particles

        !メモリの確保
        allocate(particle_position(number_of_particles,3))
        allocate(particle_velocity(number_of_particles,3))
        allocate(particle_acceleration(number_of_particles,3))
        allocate(velocity_after_collision(number_of_particles,3))
        allocate(number_density(number_of_particles))
        allocate(boundary_condition(number_of_particles))
        allocate(particle_type(number_of_particles))
        allocate(particle_pressure(number_of_particles))
        allocate(Original_layer(number_of_particles))
        allocate(source_term(number_of_particles))
        allocate(minimum_pressure(number_of_particles))

        velocity_after_collision(:,:) = 0.0
        number_density(:) = 0.0_8
        boundary_condition(:) = 0
        particle_pressure(:) = 0.0
        minimum_pressure(:) = 0.0


        !粒子の初期状態設定
        i=0
        z = 0.0-real((wallthickness+dummywallthickness),8)*particle_distance
        do while (z<=max_z+particle_distance*epsilon)
            y = 0.0-(wallthickness+dummywallthickness)*particle_distance

            do while (y<=(y_watertank+real((wallthickness+dummywallthickness),8)*particle_distance)+particle_distance*epsilon)

                x = 0.0-(wallthickness+dummywallthickness)*particle_distance
                do while (x<=x_watertank+(real((wallthickness+dummywallthickness),8)*particle_distance)+particle_distance*epsilon)

                    if (z<=z_watertank .and. (x<0 .or. y<0 .or. z<0 .or. x>x_watertank .or. y>y_watertank) ) then

                        if (x<0.0-real(wallthickness,8)*particle_distance .or. &
                        y<0.0-real(wallthickness,8)*particle_distance .or. &
                        z<0.0-real(wallthickness,8)*particle_distance .or. &
                        x>x_watertank+real(wallthickness,8)*particle_distance .or. &
                        y>y_watertank+real(wallthickness,8)*particle_distance )then
                            !ダミー粒子
                            i=i+1
                            particle_position(i,:) = (/x,y,z/)
                            particle_velocity(i,:) = 0.0
                            particle_acceleration(i,:) = 0.0
                            particle_type(i) = dummywall
                            Original_layer(i) = z
                        else
                            !壁粒子
                            i=i+1
                            particle_position(i,:) = (/x,y,z/)
                            particle_velocity(i,:) = 0.0
                            particle_acceleration(i,:) = 0.0
                            particle_type(i) = wall
                            Original_layer(i) = z
                        end if

                    else if (x>=x_watercolumn_min .and. y>=y_watercolumn_min .and. z>=z_watercolumn_min .and. x<= x_watercolumn_max .and. y<=y_watercolumn_max .and. z<= z_watercolumn_max) then
                        !流体粒子
                        i=i+1
                        particle_position(i,:) = (/x,y,z/)
                        particle_velocity(i,:) = 0.0
                        particle_acceleration(i,:) = 0.0
                        particle_type(i) = fluid
                        Original_layer(i) = z
                    end if

                    x = x + particle_distance
                end do
                y = y + particle_distance
            end do
            z = z + particle_distance
        end do



        

    end subroutine

end module initial_particle_position_velocity_particle_type

module function_module  
    implicit none
    contains 

    real(8) function weight_function(distance,Re)
        implicit none
        real(8),intent(in)::distance,Re
        if (distance<Re) then
            weight_function =(Re/distance)+(distance/Re)-2.0
        else if (distance>=Re) then
            weight_function = 0.0
        end if
        return 
    end function weight_function

    real(8) function weight_function_grad(distance,Re)
        implicit none
        real(8),intent(in)::distance,Re
        if (distance<Re) then
            weight_function_grad =(Re/distance)-(distance/Re)
        else if (distance>=Re) then
            weight_function_grad = 0.0
        end if
        return 
    end function weight_function_grad

end module function_module

module calculation_of_parameters
    use parameters_and_variables_for_simulation
    use function_module
    implicit none

    contains
    subroutine calculation_parameters
        implicit none
        call calc_n0_and_lambda()
    end subroutine calculation_parameters

    subroutine calc_n0_and_lambda
        !real(8) :: n0_for_laplacian を計算する。
        !real(8) :: lambda_0を計算する
        implicit none
        integer :: i_for_Re
        integer :: iX,iY,iZ
        integer :: iZ_max
        real(8) :: xdist,ydist,zdist,distance
        real(8) :: w,lambda,w_number_density,w_gradient


        i_for_Re = floor(Re_for_laplacian/particle_distance)+1
        n0_for_laplacian = 0.0_8
        n0_for_number_density = 0.0
        n0_for_gradient = 0.0
        lambda = 0.0_8
        
        
        if (dimention == 2) then
            !$omp parallel do reduction(+:n0_for_laplacian,n0_for_number_density,n0_for_gradient,lambda) &
            !$omp& private(iY, xdist, ydist, distance, w, w_number_density, w_gradient)
            do iX = -i_for_Re,i_for_Re  
                do iY = -i_for_Re,i_for_Re  
                    xdist = real(iX,8)*particle_distance
                    ydist = real(iY,8)*particle_distance
                    distance = sqrt(xdist**2 + ydist**2)
                    if(iX == 0 .and. iY == 0) then
                        w = 0.0
                        w_number_density = 0.0
                        w_gradient = 0.0
                    else
                        w = weight_function(distance,Re_for_laplacian) 
                        w_number_density = weight_function(distance,Re_for_number_density)
                        w_gradient = weight_function_grad(distance,Re_for_gradient)
                    end if
                    n0_for_laplacian = n0_for_laplacian + w
                    n0_for_number_density = n0_for_number_density + w_number_density
                    n0_for_gradient = n0_for_gradient + w_gradient
                    lambda = lambda + (distance**2) * w
                    
                end do
            end do
            !$omp end parallel do 
            lambda_0 = lambda/n0_for_laplacian
        end if

        if (dimention== 3) then  
            !$omp parallel do reduction(+:n0_for_laplacian,n0_for_number_density,n0_for_gradient,lambda) &
            !$omp& private(iY,iZ, xdist, ydist, zdist, distance, w, w_number_density, w_gradient)
            do iX = -i_for_Re,i_for_Re  
                do iY = -i_for_Re,i_for_Re  
                    do iZ = -i_for_Re,i_for_Re  
                        xdist = real(iX,8)*particle_distance
                        ydist = real(iY,8)*particle_distance
                        zdist = real(iZ,8)*particle_distance
                        distance = sqrt(xdist**2.0 + ydist**2.0 + zdist**2.0)
                        if(iX == 0 .and. iY == 0 .and. iZ == 0) then
                            w = 0
                            w_number_density = 0.0
                            w_gradient = 0.0
                        else
                            w = weight_function(distance,Re_for_laplacian) 
                            w_number_density = weight_function(distance,Re_for_number_density)
                            w_gradient = weight_function_grad(distance,Re_for_gradient)
                        end if
                        n0_for_laplacian = n0_for_laplacian + w
                        n0_for_number_density = n0_for_number_density + w_number_density
                        n0_for_gradient = n0_for_gradient + w_gradient
                        lambda = lambda + (distance**2) * w
                    end do
                    
                end do
            end do
            !$omp end parallel do 
            lambda_0 = lambda/n0_for_laplacian
        end if
    end subroutine calc_n0_and_lambda
end module calculation_of_parameters

module calculation_module
    use parameters_and_variables_for_simulation
    use function_module
    use mpi_module
    implicit none
    contains

    subroutine calgravity()
        !重力項による粒子にかかる加速度
        implicit none
        !内部変数
        integer :: i

        !$omp parallel do
        do i = i_start,i_finish
            if (particle_type(i)== fluid) then
                particle_acceleration(i,1) = particle_acceleration(i,1)+g_x
                particle_acceleration(i,2) = particle_acceleration(i,2)+g_y
                particle_acceleration(i,3) = particle_acceleration(i,3)+g_z
            end if
        end do
        !$omp end parallel do

    end subroutine calgravity

    subroutine calviscosity()
        implicit none
        real(8) :: a
        integer :: i,j
        real(8) :: viscosity_term(3)
        real(8) :: xdist,ydist,zdist,distance
        real(8) :: w

        a = (viscosity*2.0_8*dimention)/(n0_for_laplacian*lambda_0)


        !$omp parallel do private(j, xdist, ydist, zdist, distance, w, viscosity_term)
        loopi : do i = i_start,i_finish
            !wall,dummywallについては計算しない
            if(Particle_type(i)==wall .or. particle_type(i)==dummywall) then
                cycle loopi
            end if
            viscosity_term(:) = 0.0_8

            !i以外のjについての計算。
            loopj : do j=1,number_of_particles
                !dummywallはviscositytermの計算に利用しない
                !j=iについては和を取らない
                if(Particle_type(j)==dummywall .or. j==i) then
                    cycle loopj
                end if

                !具体的な計算
                !|xj-xi|
                xdist = particle_position(j,1)-particle_position(i,1)
                ydist = particle_position(j,2)-particle_position(i,2)
                zdist = particle_position(j,3)-particle_position(i,3)
                distance = sqrt(xdist**2.0_8 + ydist**2.0_8 + zdist**2.0_8)

                !メインのsum
                if (distance < Re_for_laplacian*1.1_8) then
                    w = weight_function(distance,Re_for_laplacian)
                    viscosity_term(1) = viscosity_term(1)+(particle_velocity(j,1)-particle_velocity(i,1))*w
                    viscosity_term(2) = viscosity_term(2)+(particle_velocity(j,2)-particle_velocity(i,2))*w
                    viscosity_term(3) = viscosity_term(3)+(particle_velocity(j,3)-particle_velocity(i,3))*w
                end if
            end do loopj

            !係数をかける
            viscosity_term(:) = a*viscosity_term(:)

            !加速度に加える
            particle_acceleration(i,1)=particle_acceleration(i,1)+viscosity_term(1)
            particle_acceleration(i,2)=particle_acceleration(i,2)+viscosity_term(2)
            particle_acceleration(i,3)=particle_acceleration(i,3)+viscosity_term(3)

        end do loopi
        !$omp end parallel do

    end subroutine calviscosity

    subroutine moveparticle()
        !加速度を受けて粒子の位置を更新する
        implicit none
        !内部変数
        integer :: i

        !$omp parallel do
        do i= i_start,i_finish
            if (particle_type(i) == fluid) then
                !速度の更新
                particle_velocity(i,1) = particle_velocity(i,1)+particle_acceleration(i,1)*time_interval
                particle_velocity(i,2) = particle_velocity(i,2)+particle_acceleration(i,2)*time_interval
                particle_velocity(i,3) = particle_velocity(i,3)+particle_acceleration(i,3)*time_interval
                !位置の更新
                particle_position(i,1) = particle_position(i,1)+particle_velocity(i,1)*time_interval
                particle_position(i,2) = particle_position(i,2)+particle_velocity(i,2)*time_interval
                particle_position(i,3) = particle_position(i,3)+particle_velocity(i,3)*time_interval
            end if

            particle_acceleration(i,:) = 0.0

        end do
        !end $omp parallel do

        call allgather_real_vector3(particle_position)
        call allgather_real_vector3(particle_velocity)
        call allgather_real_vector3(particle_acceleration)

    end subroutine moveparticle

    subroutine collision
        !粒子が衝突していると判定されたら、衝突インパルスを与え粒子の位置を更新する。
        implicit none
        !内部変数
        integer :: i,j
        real(8) :: xij,yij,zij
        real(8) :: collision_distance2
        real(8) :: distance2,distance
        real(8) :: velocity_ix, velocity_iy, velocity_iz
        real(8) :: forceDT !衝突インパルス
        real(8) :: relative_speed_negative

        collision_distance2 = collision_distance ** 2.0

        !$omp parallel do private(j, xij, yij, zij, distance2, distance, velocity_ix,velocity_iy,velocity_iz,forceDT,relative_speed_negative)
        do i= i_start,i_finish
        if (particle_type(i) == fluid) then
            velocity_ix = particle_velocity(i,1)
            velocity_iy = particle_velocity(i,2)
            velocity_iz = particle_velocity(i,3)
            velocity_after_collision(i,1) = velocity_ix
            velocity_after_collision(i,2) = velocity_iy
            velocity_after_collision(i,3) = velocity_iz

            do j=1,number_of_particles
            if (j==i) then
                cycle
            end if

                !相対距離の計算
                xij = particle_position(j,1)- particle_position(i,1)
                yij = particle_position(j,2)- particle_position(i,2)
                zij = particle_position(j,3)- particle_position(i,3)

                distance2 = (xij*xij + yij*yij + zij*zij) 

                if (distance2 < collision_distance2 .and. distance2 > 0.0) then

                    distance = sqrt(distance2)

                    !衝突向きの速度-u_ijの計算(forceDTに格納してしまって良い。)
                    relative_speed_negative &
                    = (velocity_ix-particle_velocity(j,1))*xij/distance &
                    + (velocity_iy-particle_velocity(j,2))*yij/distance &
                    + (velocity_iz-particle_velocity(j,3))*zij/distance

                    if (relative_speed_negative > 0.0) then !衝突向きの速度を持つ場合
                        forceDT = ((1.0+restitution_coefficient)/2.0)*relative_speed_negative
                        !衝突による速度の更新
                        velocity_ix = velocity_ix - forceDT*xij/distance
                        velocity_iy = velocity_iy - forceDT*yij/distance
                        velocity_iz = velocity_iz - forceDT*zij/distance

                    end if

                end if
                
            end do

            velocity_after_collision(i,1) = velocity_ix
            velocity_after_collision(i,2) = velocity_iy
            velocity_after_collision(i,3) = velocity_iz

        end if
        end do
        !end $omp parallel do

        !位置・速度の更新
        !$omp parallel do
        do i = i_start, i_finish
        if (particle_type(i) == fluid) then
            particle_position(i,1) = particle_position(i,1) + (velocity_after_collision(i,1)-particle_velocity(i,1))*time_interval
            particle_position(i,2) = particle_position(i,2) + (velocity_after_collision(i,2)-particle_velocity(i,2))*time_interval
            particle_position(i,3) = particle_position(i,3) + (velocity_after_collision(i,3)-particle_velocity(i,3))*time_interval
            particle_velocity(i,1) = velocity_after_collision(i,1)
            particle_velocity(i,2) = velocity_after_collision(i,2)
            particle_velocity(i,3) = velocity_after_collision(i,3)
        end if
        end do
        !end $omp parallel do

        call allgather_real_vector3(particle_position)
        call allgather_real_vector3(particle_velocity)

    end subroutine collision

    subroutine calnumberdensity()
        !粒子i近傍の粒子数密度を計算するルーチン
        implicit none
        !内部変数
        integer :: i,j
        real(8) :: xij,yij,zij
        real(8) :: distance,distance2

        !$omp parallel do private(j,xij,yij,zij,distance2,distance)
        do i = i_start,i_finish
            number_density(i) = 0.0
            do j =1,number_of_particles

                !粒子数密度計算
                xij = particle_position(j,1)-particle_position(i,1)
                yij = particle_position(j,2)-particle_position(i,2)
                zij = particle_position(j,3)-particle_position(i,3)

                distance2 = xij*xij+yij*yij+zij*zij
                if ( distance2 >= Re_for_number_density*Re_for_number_density .or. distance2 <= 0.0) then
                    cycle
                end if

                distance = sqrt(distance2)
                number_density(i) = number_density(i)+weight_function(distance,Re_for_number_density)

            end do
        end do
        !$omp end parallel do

        call allgather_real_vector(number_density)

    end subroutine calnumberdensity

    subroutine set_boundary_condition()
        !自由境界表面を識別するルーチン
        implicit none
        !内部変数
        integer :: i

        !$omp parallel do 
        do i = i_start,i_finish
            if (particle_type(i) == dummywall) then
                boundary_condition(i) = dummy
            else if (number_density(i)<threshold_ratio_of_number_density*n0_for_number_density) then
                boundary_condition(i) = surface_particle
            else
                boundary_condition(i) = inner_particle
            end if


        end do
        !$omp end parallel do

        call allgather_integer_vector(boundary_condition)

    end subroutine set_boundary_condition

    subroutine calPressureExplicit()
        implicit none

        integer :: i
        real(8) :: soundspeed2

        soundspeed2 = sound_speed_for_calculation* sound_speed_for_calculation
        
        !$omp parallel do
        do i = i_start,i_finish
            particle_pressure(i)=0.0
            if (boundary_condition(i) /= inner_particle) cycle
            particle_pressure(i) = soundspeed2*fluid_density*(number_density(i)-n0_for_number_density)/n0_for_number_density
        end do  
        !$omp end parallel do

    end subroutine calPressureExplicit

    subroutine remove_negative_pressure()
        implicit none
        integer :: i

        !$omp parallel do
        do i = i_start,i_finish
            if (particle_pressure(i) < 0.0_8) particle_pressure(i) = 0.0_8
        end do
        !$omp end parallel do

        call allgather_real_vector(particle_pressure)

    end subroutine remove_negative_pressure

    subroutine calpressuregradient()
        !圧力勾配から加速度を計算する。
        implicit none
        integer :: i, j
        real(8) :: gradient_x, gradient_y, gradient_z
        real(8) :: xij, yij, zij
        real(8) :: distance2, distance
        real(8) :: w, pij
        real(8) :: a
        real(8) :: re2_gradient

        a = real(dimention,8) / n0_for_gradient
        re2_gradient = Re_for_gradient * Re_for_gradient

        !$omp parallel do private(j, gradient_x, gradient_y, gradient_z, xij, yij, zij, distance2, distance, w, pij)
        do i = i_start,i_finish
            if (particle_type(i) /= fluid) then
                cycle
            end if

            gradient_x = 0.0
            gradient_y = 0.0
            gradient_z = 0.0

            do j = 1, number_of_particles
                if (j == i) then
                    cycle
                end if
                if (particle_type(j) == dummywall) then
                    cycle
                end if

                xij = particle_position(j,1) - particle_position(i,1)
                yij = particle_position(j,2) - particle_position(i,2)
                zij = particle_position(j,3) - particle_position(i,3)

                distance2 = xij*xij + yij*yij + zij*zij
                if (distance2 >= re2_gradient .or. distance2 <= 0.0) cycle

                distance = sqrt(distance2)
                w = weight_function_grad(distance, Re_for_gradient)
                pij = (particle_pressure(j) + particle_pressure(i)) / distance2

                gradient_x = gradient_x + xij * pij * w
                gradient_y = gradient_y + yij * pij * w
                gradient_z = gradient_z + zij * pij * w
            end do

            gradient_x = gradient_x * a
            gradient_y = gradient_y * a
            gradient_z = gradient_z * a

            particle_acceleration(i,1) = -gradient_x / fluid_density
            particle_acceleration(i,2) = -gradient_y / fluid_density
            particle_acceleration(i,3) = -gradient_z / fluid_density
        end do
        !$omp end parallel do

    end subroutine calpressuregradient

    subroutine moveparticleusingpressuregradient()
        !圧力勾配で得た加速度を用いて速度と位置を補正する
        implicit none
        integer :: i

        !$omp parallel do
        do i = i_start,i_finish
            if (particle_type(i) == fluid) then
                particle_velocity(i,1) = particle_velocity(i,1) + particle_acceleration(i,1) * time_interval
                particle_velocity(i,2) = particle_velocity(i,2) + particle_acceleration(i,2) * time_interval
                particle_velocity(i,3) = particle_velocity(i,3) + particle_acceleration(i,3) * time_interval

                particle_position(i,1) = particle_position(i,1) + particle_acceleration(i,1) * time_interval * time_interval
                particle_position(i,2) = particle_position(i,2) + particle_acceleration(i,2) * time_interval * time_interval
                particle_position(i,3) = particle_position(i,3) + particle_acceleration(i,3) * time_interval * time_interval
            end if

            particle_acceleration(i,:) = 0.0_8
        end do
        !$omp end parallel do

        call allgather_real_vector3(particle_position)
        call allgather_real_vector3(particle_velocity)

    end subroutine moveparticleusingpressuregradient

end module calculation_module

module output_module
    use parameters_and_variables_for_simulation
    implicit none

    contains
    subroutine writedatainvtuformat(file_number)
        !VTKファイルに粒子についての情報を出力するサブルーチン。
        !file_numberは（連番）VTKファイルの番号。
        implicit none
        character(128) :: file_name
        character(256) :: temp_char
        integer,intent(in) :: file_number
        integer :: i

        !./output_vtuファイル下の、output_(file_number).vtuに出力する。
        !現在、file_numberは６桁まで対応している。
        write(file_name,'(A,I6.6,A)') "./output_vtu/output_",file_number,".vtu"
        
        !出力部分
        open(10,file=file_name,action='write',status='replace')
        !地の文
        write(10,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
        write(10,'(A)')"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
        write(10,'(A)') "   <UnstructuredGrid>"
        write(temp_char,'(A,I0,A,I0,A)')"      <Piece NumberOfCells='",number_of_particles,"' NumberOfPoints='",number_of_particles,"'>"
        write(10,'(A)') temp_char
        !粒子の座標出力
        write(10,*) 
        write(10,'(A)') "         <Points>"
        write(10,'(A,I0,A)') "            <DataArray NumberOfComponents='3' type='Float32' Name='Particle_Position' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)particle_position(i,:)
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Points>"
        write(10,*)
        write(10,'(A)') "         <PointData>"
        !particle_typeの出力
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Int32' Name='Particle_type' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) particle_type(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !絶対速度の出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='abs_velocity' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) sqrt(particle_velocity(i,1)**2 + particle_velocity(i,2)**2 + particle_velocity(i,3)**2)
        end do
        write(10,'(A)') "            </DataArray>"
        !圧力の出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) particle_pressure(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !表面粒子か否かの出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='surface' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) boundary_condition(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !初期の液体の高さの出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='Original_layer' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) Original_layer(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !その他記述(cellなど)
        write(10,'(A)') "         </PointData>"
        write(10,*)
        write(10,'(A)') "         <Cells>"
        write(10,'(A)') "            <DataArray type='Int32' Name='connectivity' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i-1
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='Int32' Name='offsets' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='Int32' Name='types' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)1
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Cells>"
        write(10,*)
        write(10,'(A)') "      </Piece>"
        write(10,'(A)') "   </UnstructuredGrid>"
        write(10,'(A)') "</VTKFile>"
        close(10)
    end subroutine writedatainvtuformat
end module output_module

module mainloop
    use omp_lib
    use mpi_module
    use parameters_and_variables_for_simulation
    use output_module
    implicit none

    contains

    subroutine mainloopofsimulation()
        use calculation_module
        implicit none
        integer :: timestep,outputstep
        real(8) :: tstart,tfinish,tcurrent
        real(8) :: estimated_time


        !outputstep = 0
        !do while not end sim, call writedatainvtu(outputstep),
        !then, outputstep+=1

        tstart = omp_get_wtime()

        outputstep=0
        mainloop : do  timestep= 0,max_timestep
            !初回の処理
            if(timestep == 0)then
                call calnumberdensity()
                call set_boundary_condition()
                if (myrank == 0) then
                    call writedatainvtuformat(0)
                end if
                cycle mainloop
            end if
            !粘性項と重力項
            call calgravity()
            call calviscosity()
            call moveparticle()
            !衝突
            call collision()
            !圧力計算の準備
            call calnumberdensity()
            call set_boundary_condition()
            !圧力の計算
            call calPressureExplicit()
            !圧力勾配計算の前処理
            call remove_negative_pressure()
            !圧力項
            call calpressuregradient()
            call moveparticleusingpressuregradient()

            if (mod(timestep,20)==0) then
                tcurrent = omp_get_wtime()
                outputstep=outputstep+1
                estimated_time = (tcurrent-tstart)/real(timestep,8)*real(max_timestep,8)/60.0
                if (myrank==0) then
                    call writedatainvtuformat(outputstep)
                    write(*,'(f4.1,A,i4,A,f5.1,A,f5.1,A)') real(timestep,8)*100.0/real(max_timestep,8),"%, outputstep=",outputstep,", time=",(tcurrent-tstart)/60.0,"[min], estimated finish time=",estimated_time,"[min]"
                end if
            end if
        end do mainloop 

        tfinish = omp_get_wtime()
        
        if (myrank == 0) then
            write(*,'(A,f8.2,A)') 'time=', tfinish-tstart, '[s]'
        end if

    end subroutine
    

end module mainloop


program main 
    !---------modules----------!
    use mpi_f08
    use mpi_module
    use initial_particle_position_velocity_particle_type
    use calculation_of_parameters
    use mainloop
    !---------modules----------!
    implicit none
    integer :: ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    !-------calling subroutines-------!
    if (myrank == 0) then
        write(*,*) "-----Starting initialization-----"
    end if
    call water_tank_and_water_column_2d(real(2.0,8),real(1.2,8),real(0.5,8),real(1.2 ,8),3,2)
    if (myrank ==0 ) then
        write(*,*) "-----Done initialization-----"
        write(*,*)
    end if
    call calculation_parameters
    call setup_mpi_decomposition(number_of_particles)
    call mainloopofsimulation
    !-------calling subroutines-------!

    call MPI_Finalize(ierr)


end program