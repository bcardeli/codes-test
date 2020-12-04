module light_competition
    use ISO_FORTRAN_ENV, only: REAL32, REAL64, REAL128
    use allometry

    implicit none

    type :: layer_array
        real (REAL64) :: sum_height
        integer :: num_height !!corresponds to the number of pls
        real (REAL64) :: mean_height
        real (REAL64) :: layer_height
        real (REAL64) :: sum_LAI !LAI sum in a layer
        real (REAL64) :: mean_LAI !mean LAI in a layer
        real (REAL64) :: beers_law !layer's light extinction
        real (REAL64) :: linc !layer's light incidence
        real (REAL64) :: lused !layer's light used (relates to light extinction - Beers Law)
        real (REAL64) :: lavai !light availability
    end type layer_array

    integer, parameter :: npls = 14 !number of PLS to test the logic and dynamic.
    
    
    integer :: num_layer
    real (REAL64) :: max_height
    real (REAL64) :: layer_size
    real (REAL64) :: incidence_rad !Incidence radiation (relates do APAR) in J/m2/s
    real (REAL64) :: watt_rs = 210 !shortwave radiation in watts/m2
    real (REAL64) :: short_rad !shortwave radiation in joules/s

    integer::i,j

    integer :: last_with_pls

    type(layer_array), allocatable :: layer(:)

    ! Layer's dynamics functions

    max_height = maxval(height)
    !print*, 'max_height',max_height
    
    num_layer = nint(max_height/5)
    !print*, 'num_layer',num_layer

    allocate(layer(1:num_layer))
    
    last_with_pls=num_layer

    layer_size = max_height/num_layer
    !print*, 'layer_size', layer_size

    layer(i)%layer_height=0

    do i=1,num_layer
        layer(i)%layer_height=layer_size*i
        !print*, 'layer_height',layer(i)%layer_height, i
    enddo


end module light_competition