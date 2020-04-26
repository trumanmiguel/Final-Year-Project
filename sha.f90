! Goal is to calculate the spin Hall angle given the spin Hall conductivity and the charge conductivity
! Change nfermi, seedname, MuNumPoints, TempNumPoints, shc_beta, fermi_energy

Program sha
    implicit none
    integer :: nfermi = 3, i, MuIdx, TempIdx
    real :: fermi_energy = 8.7816
    real, allocatable, dimension(:,:) :: spinHA
    character(len=50) :: seedname = "WTe2"

    ! Parameters for spin Hall conductivity (SHC): units of (hbar/e) S/cm
    real, allocatable, dimension(:) :: fermi_energy_list(:), shc_fermi(:)
    integer :: shc_beta = 2, elec_field = 0

    ! Parameters for electrical conductivity (EC): units of 1/(Ohm*m) = S/m
    integer :: MuNumPoints = 12, TempNumPoints = 7
    real, allocatable, dimension(:,:,:) :: ElCond ! (TE (tensor elements; W90 developers call TE as coordinate) ,Temp, Mu)
    ! 1 = xx, 2 = xy, 3 = yy, 4 = xz, 5 = yz, 6 = zz
    real, allocatable, dimension(:) :: TempArray, MuArray

    ! Allocation for SHC related arrays
    allocate (fermi_energy_list(nfermi))
    allocate (shc_fermi(nfermi))

    ! Allocation for EC related arrays
    allocate (ElCond(6,TempNumPoints,MuNumPoints))
    allocate (MuArray(MuNumPoints)) ! In general, MuArray should be the same size as shc_fermi
    allocate (TempArray(TempNumPoints))

    ! Allocation for SHA related arrays
    allocate (spinHA(nfermi,TempNumPoints))

    ! Reading of SHC data
    open(1, file = trim(seedname)//'-shc-fermiscan.dat', status = 'old')

    read(1,*) ! "Read" the first header line by doing nothing

    ! Read only the Fermi energy and the SHC
    do i = 1,nfermi
        read(1,101) fermi_energy_list(i), shc_fermi(i)
    end do

    101 format (8x,f9.6,2x,f15.3)

    close(1)

    ! Reading of EC data
    open(2, file = trim(seedname)//'_elcond.dat', status = 'old')

    do i = 1,3
        read(2,*) ! "Read" the first 3 header lines by doing nothing
    end do

    ! Read the chemical potential (eV), temperature (K) and the EC tensor elements
    do MuIdx = 1,MuNumPoints
        do TempIdx = 1,TempNumPoints
            read(2,*) MuArray(MuIdx), TempArray(TempIdx), ElCond(:,TempIdx,MuIdx)
        end do
    end do

    close(2)

    ! Calculation of the spin Hall angle (SHA)
    ! Spin Hall angle equation is taken from Nature Materials volume 19, pages 292–298(2020) under the methods section
    open(3, file = trim(seedname)//'-sha-fermiscan.dat', status = 'replace')
    write (3,103) 'Fermi energy(eV)', 'Temp(K)', 'SHA'

    ! Choosing which electrical conductivity tensor element to divide by in the spin Hall angle equation
    select case (shc_beta)
        case (1)
            elec_field = 1
        case (2)
            elec_field = 3
        case (3)
            elec_field = 6
    end select

    ! Calculating the spin Hall angle matrix. Rows: Fermi energy; Columns: Temperature
    do MuIdx = 1,MuNumPoints
        do TempIdx = 1,TempNumPoints
            do i = 1,nfermi
                if (fermi_energy_list(i) == MuArray(MuIdx)) then
                    spinHA(i,TempIdx) = 2*((shc_fermi(i)/0.01)/ElCond(elec_field,TempIdx,MuIdx))
                    write(3,*), fermi_energy_list(i), TempArray(TempIdx), spinHA(i,TempIdx)
                end if
            end do
        end do
    end do

    103 format (a,5x,a,12x,a)

    close(3)

End Program

