module tt_bbeevv

    use kinds

    implicit none

    real(kind=default), public :: sgg_tt_bbeevv
    real(kind=default), private :: gg_tt_bbeevv
    real(kind=default), public :: sqq_tt_bbeevv
    real(kind=default), private :: qq_tt_bbeevv
    real(kind=default), public :: sqq_tt_bbeevv_ew
    real(kind=default), private :: qq_tt_bbeevv_ew

contains

function sgg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8, channel)

    ! returns amplitude squared summed/avg over colors and helicities
    ! for the point in phase space p1, p2, p3, p4, ...
    ! for process: g g -> t t~ -> b b~ ta+ ta- vt vt~

    implicit none

    ! constants
    integer, parameter :: nexternal = 8, ncomb = 256

    ! arguments
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
    integer, intent(in), optional :: channel

    ! local variables
    integer :: i, j
    integer :: nhel(nexternal, ncomb), ntry
    real(kind=default) :: t
    integer :: ihel
    logical :: goodhel(ncomb)
    data goodhel /ncomb * .false./
    data ntry /0/

    ! possible helicity combinations
    data (nhel(ihel,   1), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel,   2), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel,   3), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel,   4), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel,   5), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel,   6), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel,   7), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel,   8), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel,   9), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  10), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  11), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  12), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  13), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  14), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  15), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  16), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  17), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  18), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  19), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  20), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  21), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  22), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  23), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  24), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  25), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  26), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  27), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  28), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  29), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  30), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  31), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  32), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  33), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  34), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  35), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel,  36), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel,  37), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel,  38), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel,  39), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel,  40), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel,  41), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  42), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  43), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  44), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  45), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  46), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  47), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  48), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  49), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  50), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  51), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  52), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  53), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  54), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  55), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  56), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  57), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  58), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  59), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  60), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  61), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  62), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  63), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  64), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  65), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  66), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  67), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel,  68), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel,  69), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel,  70), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel,  71), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel,  72), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel,  73), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  74), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  75), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  76), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  77), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  78), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  79), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  80), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  81), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  82), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  83), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  84), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  85), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  86), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  87), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  88), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  89), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  90), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  91), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  92), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  93), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  94), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  95), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  96), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  97), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  98), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  99), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 100), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 101), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 102), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 103), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 104), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 105), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 106), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 107), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 108), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 109), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 110), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 111), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 112), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 113), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 114), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 115), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 116), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 117), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 118), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 119), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 120), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 121), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 122), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 123), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 124), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 125), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 126), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 127), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 128), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 129), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 130), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 131), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 132), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 133), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 134), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 135), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 136), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 137), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 138), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 139), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 140), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 141), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 142), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 143), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 144), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 145), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 146), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 147), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 148), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 149), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 150), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 151), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 152), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 153), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 154), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 155), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 156), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 157), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 158), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 159), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 160), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 161), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 162), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 163), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 164), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 165), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 166), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 167), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 168), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 169), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 170), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 171), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 172), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 173), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 174), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 175), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 176), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 177), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 178), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 179), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 180), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 181), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 182), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 183), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 184), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 185), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 186), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 187), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 188), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 189), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 190), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 191), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 192), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 193), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 194), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 195), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 196), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 197), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 198), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 199), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 200), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 201), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 202), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 203), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 204), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 205), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 206), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 207), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 208), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 209), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 210), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 211), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 212), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 213), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 214), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 215), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 216), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 217), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 218), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 219), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 220), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 221), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 222), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 223), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 224), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 225), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 226), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 227), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 228), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 229), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 230), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 231), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 232), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 233), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 234), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 235), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 236), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 237), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 238), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 239), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 240), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 241), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 242), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 243), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 244), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 245), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 246), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 247), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 248), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 249), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 250), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 251), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 252), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 253), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 254), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 255), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 256), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1,  1/

    sgg_tt_bbeevv = 0.d0
    ntry = ntry + 1
    do ihel = 1, ncomb
        if (goodhel(ihel) .or. ntry < 10) then
            t = gg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8, nhel(1, ihel), channel)
            sgg_tt_bbeevv = sgg_tt_bbeevv + t
            if (t > 0.d0 .and. .not. goodhel(ihel)) then
                goodhel(ihel)= .true.
            endif
        end if
    enddo
    sgg_tt_bbeevv = sgg_tt_bbeevv / 4.d0
end function sgg_tt_bbeevv


function gg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8, nhel, channel)

    ! returns amplitude squared summed/avg over colors
    ! for the point in phase space p1, p2, p3, p4, ...
    ! and helicity nhel(1), nhel(2), ...
    ! for process : g g  -> t t~ -> b b~ ta+ ta- vt vt~

    use helas
    use modelling

    implicit none

    ! constants
    integer, parameter :: ngraphs = 3, neigen = 2, nexternal = 8
    real(kind=default), parameter :: zero = 0.d0

    ! arguments
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
    integer :: nhel(nexternal)
    integer, intent(in), optional :: channel

    ! local variables
    integer :: i, j
    real(kind=default) :: eigen_val(neigen), eigen_vec(ngraphs, neigen)
    complex(kind=complex) :: atmp, atmp2
    complex(kind=complex) :: amp(ngraphs)
    complex(kind=complex) :: w1(6), w2(6), w3(6), w4(6), w5(6)
    complex(kind=complex) :: w6(6), w7(6), w8(6), w9(6), w10(6)
    complex(kind=complex) :: w11(6), w12(6), w13(6), w14(6), w15(6)
    real(kind=default) :: atot2, asum2

    ! initialise amplitudes
    do i = 1, ngraphs
        amp(i) = 0.d0
    enddo

    ! color data
    eigen_val(1) = 7.2916666666666588e-02
    eigen_vec(1, 1) = 7.0710678118654768e-01
    eigen_vec(2, 1) = 7.0710678118654735e-01
    eigen_vec(3, 1) = 0.0000000000000000e+00
    eigen_val(2) = 2.8125000000000000e-01
    eigen_vec(1, 2) = -4.0824829046386313e-01
    eigen_vec(2, 2) = 4.0824829046386285e-01
    eigen_vec(3, 2) = 8.1649658092772603e-01

    ! wavefunctions
    call vxxxxx(p1, zero, nhel(1), -1, w1)
    call vxxxxx(p2, zero, nhel(2), -1, w2)
    call oxxxxx(p3, fmass(12), nhel(3), 1, w3)
    call ixxxxx(p4, fmass(12), nhel(4), -1, w4)
    call ixxxxx(p5, zero, nhel(5), -1, w5)
    call oxxxxx(p6, zero, nhel(6), 1, w6)
    call oxxxxx(p7, zero, nhel(7), 1, w7)
    call ixxxxx(p8, zero, nhel(8), -1, w8)

    ! currents
    call jioxxx(w5, w7, gwf, wmass, wwidth, w9)
    call jioxxx(w8, w6, gwf, wmass, wwidth, w10)
    call fvoxxx(w3, w9, gwf, fmass(11), fwidth(11), w11)
    call fvixxx(w4, w10, gwf, fmass(11), fwidth(11), w12)
    call fvoxxx(w11, w2, gg, fmass(11), fwidth(11), w13)
    call fvoxxx(w11, w1, gg, fmass(11), fwidth(11), w14)
    call jggxxx(w1, w2, g, w15)

    ! amplitudes
    call iovxxx(w12, w13, w1, gg, amp(1))
    call iovxxx(w12, w14, w2, gg, amp(2))
    call iovxxx(w12, w11, w15, gg, amp(3))

    ! full square matrix element
    atot2 = 0.d0
    do i = 1, neigen
        atmp = (0.d0, 0.d0)
        do j = 1, ngraphs
            atmp = atmp + eigen_vec(j,i) * amp(j)
        enddo
        atot2 = atot2 + atmp * eigen_val(i) * conjg(atmp)
    enddo

    if (present(channel) .and. channel > 0) then

        ! sum of individual square amplitudes
        asum2 = 0.d0
        do i = 1, neigen
            atmp = (0.d0, 0.d0)
            do j = 1, ngraphs
                atmp = eigen_vec(j,i) * amp(j)
                asum2 = asum2 + atmp * eigen_val(i) * conjg(atmp)
            enddo
        enddo

        ! individual square amplitude for this channel
        if (asum2 > 0.d0) then
            atmp2 = 0.d0
            do i = 1, neigen
                atmp = eigen_vec(channel, i) * amp(channel)
                atmp2 = atmp2 + atmp * eigen_val(i) * conjg(atmp) * atot2 / asum2
            end do
            gg_tt_bbeevv = atmp2
        else
            gg_tt_bbeevv = 0
        end if
    else
        gg_tt_bbeevv = atot2
    end if

end function gg_tt_bbeevv

function sqq_tt_bbeevv(iq, p1, p2, p3, p4, p5, p6, p7, p8)

    ! function generated by madgraph
    ! returns amplitude squared summed/avg over colors and helicities
    ! for the point in phase space p1, p2, p3, p4, ...
    ! for process: q q~ -> b b~ e+ e- v v~

    implicit none

    ! constants
    integer, parameter :: nexternal = 8, ncomb = 256

    ! arguments
    integer :: iq
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

    ! local variables
    integer :: nhel(nexternal, ncomb), ntry
    real(kind=default) :: t
    integer :: ihel
    logical :: goodhel(ncomb)
    data goodhel /ncomb*.false./
    data ntry /0/

    ! possible helicity amplitudes
    data (nhel(ihel,   1), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel,   2), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel,   3), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel,   4), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel,   5), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel,   6), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel,   7), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel,   8), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel,   9), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  10), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  11), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  12), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  13), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  14), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  15), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  16), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  17), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  18), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  19), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  20), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  21), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  22), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  23), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  24), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  25), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  26), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  27), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  28), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  29), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  30), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  31), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  32), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  33), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  34), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  35), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel,  36), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel,  37), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel,  38), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel,  39), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel,  40), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel,  41), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  42), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  43), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  44), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  45), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  46), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  47), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  48), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  49), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  50), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  51), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  52), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  53), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  54), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  55), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  56), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  57), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  58), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  59), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  60), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  61), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  62), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  63), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  64), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  65), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  66), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  67), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel,  68), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel,  69), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel,  70), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel,  71), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel,  72), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel,  73), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel,  74), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel,  75), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel,  76), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel,  77), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel,  78), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel,  79), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel,  80), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel,  81), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel,  82), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel,  83), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel,  84), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel,  85), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel,  86), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel,  87), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel,  88), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel,  89), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel,  90), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel,  91), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel,  92), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel,  93), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel,  94), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel,  95), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel,  96), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel,  97), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel,  98), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel,  99), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 100), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 101), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 102), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 103), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 104), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 105), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 106), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 107), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 108), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 109), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 110), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 111), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 112), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 113), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 114), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 115), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 116), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 117), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 118), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 119), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 120), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 121), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 122), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 123), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 124), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 125), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 126), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 127), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 128), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 129), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 130), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 131), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 132), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 133), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 134), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 135), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 136), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 137), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 138), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 139), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 140), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 141), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 142), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 143), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 144), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 145), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 146), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 147), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 148), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 149), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 150), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 151), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 152), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 153), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 154), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 155), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 156), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 157), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 158), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 159), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 160), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 161), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 162), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 163), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 164), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 165), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 166), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 167), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 168), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 169), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 170), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 171), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 172), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 173), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 174), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 175), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 176), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 177), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 178), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 179), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 180), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 181), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 182), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 183), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 184), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 185), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 186), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 187), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 188), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 189), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 190), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 191), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 192), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 193), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 194), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 195), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 196), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 197), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 198), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 199), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 200), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 201), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 202), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 203), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 204), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 205), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 206), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 207), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 208), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 209), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 210), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 211), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 212), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 213), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 214), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 215), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 216), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 217), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 218), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 219), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 220), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 221), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 222), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 223), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 224), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1,  1/
    data (nhel(ihel, 225), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1, -1/
    data (nhel(ihel, 226), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1,  1/
    data (nhel(ihel, 227), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1, -1/
    data (nhel(ihel, 228), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1,  1/
    data (nhel(ihel, 229), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1, -1/
    data (nhel(ihel, 230), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1,  1/
    data (nhel(ihel, 231), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1, -1/
    data (nhel(ihel, 232), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1,  1/
    data (nhel(ihel, 233), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1, -1/
    data (nhel(ihel, 234), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1,  1/
    data (nhel(ihel, 235), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1, -1/
    data (nhel(ihel, 236), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1,  1/
    data (nhel(ihel, 237), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1, -1/
    data (nhel(ihel, 238), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1,  1/
    data (nhel(ihel, 239), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1, -1/
    data (nhel(ihel, 240), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1,  1/
    data (nhel(ihel, 241), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1, -1/
    data (nhel(ihel, 242), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1,  1/
    data (nhel(ihel, 243), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1, -1/
    data (nhel(ihel, 244), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1,  1/
    data (nhel(ihel, 245), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1, -1/
    data (nhel(ihel, 246), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1,  1/
    data (nhel(ihel, 247), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1, -1/
    data (nhel(ihel, 248), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1,  1/
    data (nhel(ihel, 249), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1, -1/
    data (nhel(ihel, 250), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1,  1/
    data (nhel(ihel, 251), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1, -1/
    data (nhel(ihel, 252), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1,  1/
    data (nhel(ihel, 253), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1, -1/
    data (nhel(ihel, 254), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1,  1/
    data (nhel(ihel, 255), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1, -1/
    data (nhel(ihel, 256), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1,  1/
    sqq_tt_bbeevv = 0.d0
    ntry = ntry + 1
    do ihel = 1, ncomb
        if (goodhel(ihel) .or. ntry < 10) then
            t = qq_tt_bbeevv(iq, p1, p2, p3, p4, p5, p6, p7, p8, nhel(1,ihel))
            sqq_tt_bbeevv = sqq_tt_bbeevv + t
            if (t > 0d0 .and. .not. goodhel(ihel)) then
                goodhel(ihel)= .true.
            endif
        end if
    enddo
    sqq_tt_bbeevv = sqq_tt_bbeevv / 4.d0
end function sqq_tt_bbeevv


function qq_tt_bbeevv(iq, p1, p2, p3, p4, p5, p6, p7, p8, nhel)

    ! function generated by madgraph
    ! returns amplitude squared summed/avg over colors
    ! for the point in phase space p1, p2, p3, p4, ...
    ! and helicity nhel(1), nhel(2), ...
    ! for process : q q~ t t~ -> b b~ ta+ ta- vt vt~

    use helas
    use modelling

    implicit none

    ! constants
    integer :: iq
    integer, parameter :: ngraphs = 1, neigen = 1, nexternal = 8
    real(kind=default), parameter :: zero = 0d0

    ! arguments
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
    integer :: nhel(nexternal)

    ! local variables
    integer :: i, j
    real(kind=default) :: eigen_val(neigen), eigen_vec(ngraphs,neigen)
    complex(kind=complex) ztemp
    complex(kind=complex) amp(ngraphs)
    complex(kind=complex) w1(6), w2(6), w3(6), w4(6), w5(6)
    complex(kind=complex) w6(6), w7(6), w8(6), w9(6), w10(6)
    complex(kind=complex) w11(6), w12(6), w13(6)

    ! color data
    data eigen_val(1  ) /2.2222222222222221e-01/
    data eigen_vec(1, 1) /-1.0000000000000000e+00/

    ! wavefunctions
    call ixxxxx(p1, fmass(iq), nhel(1), 1, w1)
    call oxxxxx(p2, fmass(iq), nhel(2), -1, w2)
    call oxxxxx(p3, fmass(12), nhel(3), 1, w3)
    call ixxxxx(p4, fmass(12), nhel(4), -1, w4)
    call ixxxxx(p5, zero, nhel(5), -1, w5)
    call oxxxxx(p6, zero, nhel(6), 1, w6)
    call oxxxxx(p7, zero, nhel(7), 1, w7)
    call ixxxxx(p8, zero, nhel(8), -1, w8)

    ! currents
    call jioxxx(w1, w2, gg, zero, zero, w9)
    call jioxxx(w5, w7, gwf, wmass, wwidth, w10)
    call jioxxx(w8, w6, gwf, wmass, wwidth, w11)
    call fvoxxx(w3, w10, gwf, fmass(11), fwidth(11), w12)
    call fvixxx(w4, w11, gwf, fmass(11), fwidth(11), w13)

    ! amplitude
    call iovxxx(w13, w12, w9, gg, amp(1))

    qq_tt_bbeevv = 0.d0
    do i = 1, neigen
        ztemp = (0.d0, 0.d0)
        do j = 1, ngraphs
            ztemp = ztemp + eigen_vec(j,i) * amp(j)
        enddo
        qq_tt_bbeevv = qq_tt_bbeevv + ztemp*eigen_val(i) * conjg(ztemp)
    enddo
    ! call gaugecheck(amp, ztemp, eigen_vec, eigen_val, ngraphs, neigen)
end function qq_tt_bbeevv

function sqq_tt_bbeevv_ew(iq, jf, p1, p2, p3, p4, p5, p6, p7, p8, channel)

  ! file originally generated by madgraph
  ! returns amplitude squared summed/avg over colors and helicities
  ! for the point in phase space p1, p2, p3, p4, p5, p6, p7, p8
  ! for process: q q -> A, Z, Z' -> b b~ e+ e- v v~

  use kinds

  implicit none

  ! constants
  integer, parameter :: nexternal = 8, ncomb = 256

  ! arguments
  integer :: iq, jf
  real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
  integer, intent(in), optional :: channel

  ! local variables
  integer :: nhel(nexternal, ncomb), ntry
  real(kind=default) :: t
  integer :: ihel
  logical :: goodhel(ncomb)
  data goodhel /ncomb*.false./
  data ntry /0/

  ! possible helicity combinations
  data (nhel(ihel,  1), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1, -1/
  data (nhel(ihel,  2), ihel = 1, 8) /-1, -1, -1, -1, -1, -1, -1,  1/
  data (nhel(ihel,  3), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1, -1/
  data (nhel(ihel,  4), ihel = 1, 8) /-1, -1, -1, -1, -1, -1,  1,  1/
  data (nhel(ihel,  5), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1, -1/
  data (nhel(ihel,  6), ihel = 1, 8) /-1, -1, -1, -1, -1,  1, -1,  1/
  data (nhel(ihel,  7), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1, -1/
  data (nhel(ihel,  8), ihel = 1, 8) /-1, -1, -1, -1, -1,  1,  1,  1/
  data (nhel(ihel,  9), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1, -1/
  data (nhel(ihel, 10), ihel = 1, 8) /-1, -1, -1, -1,  1, -1, -1,  1/
  data (nhel(ihel, 11), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1, -1/
  data (nhel(ihel, 12), ihel = 1, 8) /-1, -1, -1, -1,  1, -1,  1,  1/
  data (nhel(ihel, 13), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1, -1/
  data (nhel(ihel, 14), ihel = 1, 8) /-1, -1, -1, -1,  1,  1, -1,  1/
  data (nhel(ihel, 15), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1, -1/
  data (nhel(ihel, 16), ihel = 1, 8) /-1, -1, -1, -1,  1,  1,  1,  1/
  data (nhel(ihel, 17), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1, -1/
  data (nhel(ihel, 18), ihel = 1, 8) /-1, -1, -1,  1, -1, -1, -1,  1/
  data (nhel(ihel, 19), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1, -1/
  data (nhel(ihel, 20), ihel = 1, 8) /-1, -1, -1,  1, -1, -1,  1,  1/
  data (nhel(ihel, 21), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1, -1/
  data (nhel(ihel, 22), ihel = 1, 8) /-1, -1, -1,  1, -1,  1, -1,  1/
  data (nhel(ihel, 23), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1, -1/
  data (nhel(ihel, 24), ihel = 1, 8) /-1, -1, -1,  1, -1,  1,  1,  1/
  data (nhel(ihel, 25), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1, -1/
  data (nhel(ihel, 26), ihel = 1, 8) /-1, -1, -1,  1,  1, -1, -1,  1/
  data (nhel(ihel, 27), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1, -1/
  data (nhel(ihel, 28), ihel = 1, 8) /-1, -1, -1,  1,  1, -1,  1,  1/
  data (nhel(ihel, 29), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1, -1/
  data (nhel(ihel, 30), ihel = 1, 8) /-1, -1, -1,  1,  1,  1, -1,  1/
  data (nhel(ihel, 31), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1, -1/
  data (nhel(ihel, 32), ihel = 1, 8) /-1, -1, -1,  1,  1,  1,  1,  1/
  data (nhel(ihel, 33), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1, -1/
  data (nhel(ihel, 34), ihel = 1, 8) /-1, -1,  1, -1, -1, -1, -1,  1/
  data (nhel(ihel, 35), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1, -1/
  data (nhel(ihel, 36), ihel = 1, 8) /-1, -1,  1, -1, -1, -1,  1,  1/
  data (nhel(ihel, 37), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1, -1/
  data (nhel(ihel, 38), ihel = 1, 8) /-1, -1,  1, -1, -1,  1, -1,  1/
  data (nhel(ihel, 39), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1, -1/
  data (nhel(ihel, 40), ihel = 1, 8) /-1, -1,  1, -1, -1,  1,  1,  1/
  data (nhel(ihel, 41), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1, -1/
  data (nhel(ihel, 42), ihel = 1, 8) /-1, -1,  1, -1,  1, -1, -1,  1/
  data (nhel(ihel, 43), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1, -1/
  data (nhel(ihel, 44), ihel = 1, 8) /-1, -1,  1, -1,  1, -1,  1,  1/
  data (nhel(ihel, 45), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1, -1/
  data (nhel(ihel, 46), ihel = 1, 8) /-1, -1,  1, -1,  1,  1, -1,  1/
  data (nhel(ihel, 47), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1, -1/
  data (nhel(ihel, 48), ihel = 1, 8) /-1, -1,  1, -1,  1,  1,  1,  1/
  data (nhel(ihel, 49), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1, -1/
  data (nhel(ihel, 50), ihel = 1, 8) /-1, -1,  1,  1, -1, -1, -1,  1/
  data (nhel(ihel, 51), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1, -1/
  data (nhel(ihel, 52), ihel = 1, 8) /-1, -1,  1,  1, -1, -1,  1,  1/
  data (nhel(ihel, 53), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1, -1/
  data (nhel(ihel, 54), ihel = 1, 8) /-1, -1,  1,  1, -1,  1, -1,  1/
  data (nhel(ihel, 55), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1, -1/
  data (nhel(ihel, 56), ihel = 1, 8) /-1, -1,  1,  1, -1,  1,  1,  1/
  data (nhel(ihel, 57), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1, -1/
  data (nhel(ihel, 58), ihel = 1, 8) /-1, -1,  1,  1,  1, -1, -1,  1/
  data (nhel(ihel, 59), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1, -1/
  data (nhel(ihel, 60), ihel = 1, 8) /-1, -1,  1,  1,  1, -1,  1,  1/
  data (nhel(ihel, 61), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1, -1/
  data (nhel(ihel, 62), ihel = 1, 8) /-1, -1,  1,  1,  1,  1, -1,  1/
  data (nhel(ihel, 63), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1, -1/
  data (nhel(ihel, 64), ihel = 1, 8) /-1, -1,  1,  1,  1,  1,  1,  1/
  data (nhel(ihel, 65), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1, -1/
  data (nhel(ihel, 66), ihel = 1, 8) /-1,  1, -1, -1, -1, -1, -1,  1/
  data (nhel(ihel, 67), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1, -1/
  data (nhel(ihel, 68), ihel = 1, 8) /-1,  1, -1, -1, -1, -1,  1,  1/
  data (nhel(ihel, 69), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1, -1/
  data (nhel(ihel, 70), ihel = 1, 8) /-1,  1, -1, -1, -1,  1, -1,  1/
  data (nhel(ihel, 71), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1, -1/
  data (nhel(ihel, 72), ihel = 1, 8) /-1,  1, -1, -1, -1,  1,  1,  1/
  data (nhel(ihel, 73), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1, -1/
  data (nhel(ihel, 74), ihel = 1, 8) /-1,  1, -1, -1,  1, -1, -1,  1/
  data (nhel(ihel, 75), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1, -1/
  data (nhel(ihel, 76), ihel = 1, 8) /-1,  1, -1, -1,  1, -1,  1,  1/
  data (nhel(ihel, 77), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1, -1/
  data (nhel(ihel, 78), ihel = 1, 8) /-1,  1, -1, -1,  1,  1, -1,  1/
  data (nhel(ihel, 79), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1, -1/
  data (nhel(ihel, 80), ihel = 1, 8) /-1,  1, -1, -1,  1,  1,  1,  1/
  data (nhel(ihel, 81), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1, -1/
  data (nhel(ihel, 82), ihel = 1, 8) /-1,  1, -1,  1, -1, -1, -1,  1/
  data (nhel(ihel, 83), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1, -1/
  data (nhel(ihel, 84), ihel = 1, 8) /-1,  1, -1,  1, -1, -1,  1,  1/
  data (nhel(ihel, 85), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1, -1/
  data (nhel(ihel, 86), ihel = 1, 8) /-1,  1, -1,  1, -1,  1, -1,  1/
  data (nhel(ihel, 87), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1, -1/
  data (nhel(ihel, 88), ihel = 1, 8) /-1,  1, -1,  1, -1,  1,  1,  1/
  data (nhel(ihel, 89), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1, -1/
  data (nhel(ihel, 90), ihel = 1, 8) /-1,  1, -1,  1,  1, -1, -1,  1/
  data (nhel(ihel, 91), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1, -1/
  data (nhel(ihel, 92), ihel = 1, 8) /-1,  1, -1,  1,  1, -1,  1,  1/
  data (nhel(ihel, 93), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1, -1/
  data (nhel(ihel, 94), ihel = 1, 8) /-1,  1, -1,  1,  1,  1, -1,  1/
  data (nhel(ihel, 95), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1, -1/
  data (nhel(ihel, 96), ihel = 1, 8) /-1,  1, -1,  1,  1,  1,  1,  1/
  data (nhel(ihel, 97), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1, -1/
  data (nhel(ihel, 98), ihel = 1, 8) /-1,  1,  1, -1, -1, -1, -1,  1/
  data (nhel(ihel, 99), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1, -1/
  data (nhel(ihel,100), ihel = 1, 8) /-1,  1,  1, -1, -1, -1,  1,  1/
  data (nhel(ihel,101), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1, -1/
  data (nhel(ihel,102), ihel = 1, 8) /-1,  1,  1, -1, -1,  1, -1,  1/
  data (nhel(ihel,103), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1, -1/
  data (nhel(ihel,104), ihel = 1, 8) /-1,  1,  1, -1, -1,  1,  1,  1/
  data (nhel(ihel,105), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1, -1/
  data (nhel(ihel,106), ihel = 1, 8) /-1,  1,  1, -1,  1, -1, -1,  1/
  data (nhel(ihel,107), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1, -1/
  data (nhel(ihel,108), ihel = 1, 8) /-1,  1,  1, -1,  1, -1,  1,  1/
  data (nhel(ihel,109), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1, -1/
  data (nhel(ihel,110), ihel = 1, 8) /-1,  1,  1, -1,  1,  1, -1,  1/
  data (nhel(ihel,111), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1, -1/
  data (nhel(ihel,112), ihel = 1, 8) /-1,  1,  1, -1,  1,  1,  1,  1/
  data (nhel(ihel,113), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1, -1/
  data (nhel(ihel,114), ihel = 1, 8) /-1,  1,  1,  1, -1, -1, -1,  1/
  data (nhel(ihel,115), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1, -1/
  data (nhel(ihel,116), ihel = 1, 8) /-1,  1,  1,  1, -1, -1,  1,  1/
  data (nhel(ihel,117), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1, -1/
  data (nhel(ihel,118), ihel = 1, 8) /-1,  1,  1,  1, -1,  1, -1,  1/
  data (nhel(ihel,119), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1, -1/
  data (nhel(ihel,120), ihel = 1, 8) /-1,  1,  1,  1, -1,  1,  1,  1/
  data (nhel(ihel,121), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1, -1/
  data (nhel(ihel,122), ihel = 1, 8) /-1,  1,  1,  1,  1, -1, -1,  1/
  data (nhel(ihel,123), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1, -1/
  data (nhel(ihel,124), ihel = 1, 8) /-1,  1,  1,  1,  1, -1,  1,  1/
  data (nhel(ihel,125), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1, -1/
  data (nhel(ihel,126), ihel = 1, 8) /-1,  1,  1,  1,  1,  1, -1,  1/
  data (nhel(ihel,127), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1, -1/
  data (nhel(ihel,128), ihel = 1, 8) /-1,  1,  1,  1,  1,  1,  1,  1/
  data (nhel(ihel,129), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1, -1/
  data (nhel(ihel,130), ihel = 1, 8) / 1, -1, -1, -1, -1, -1, -1,  1/
  data (nhel(ihel,131), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1, -1/
  data (nhel(ihel,132), ihel = 1, 8) / 1, -1, -1, -1, -1, -1,  1,  1/
  data (nhel(ihel,133), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1, -1/
  data (nhel(ihel,134), ihel = 1, 8) / 1, -1, -1, -1, -1,  1, -1,  1/
  data (nhel(ihel,135), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1, -1/
  data (nhel(ihel,136), ihel = 1, 8) / 1, -1, -1, -1, -1,  1,  1,  1/
  data (nhel(ihel,137), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1, -1/
  data (nhel(ihel,138), ihel = 1, 8) / 1, -1, -1, -1,  1, -1, -1,  1/
  data (nhel(ihel,139), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1, -1/
  data (nhel(ihel,140), ihel = 1, 8) / 1, -1, -1, -1,  1, -1,  1,  1/
  data (nhel(ihel,141), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1, -1/
  data (nhel(ihel,142), ihel = 1, 8) / 1, -1, -1, -1,  1,  1, -1,  1/
  data (nhel(ihel,143), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1, -1/
  data (nhel(ihel,144), ihel = 1, 8) / 1, -1, -1, -1,  1,  1,  1,  1/
  data (nhel(ihel,145), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1, -1/
  data (nhel(ihel,146), ihel = 1, 8) / 1, -1, -1,  1, -1, -1, -1,  1/
  data (nhel(ihel,147), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1, -1/
  data (nhel(ihel,148), ihel = 1, 8) / 1, -1, -1,  1, -1, -1,  1,  1/
  data (nhel(ihel,149), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1, -1/
  data (nhel(ihel,150), ihel = 1, 8) / 1, -1, -1,  1, -1,  1, -1,  1/
  data (nhel(ihel,151), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1, -1/
  data (nhel(ihel,152), ihel = 1, 8) / 1, -1, -1,  1, -1,  1,  1,  1/
  data (nhel(ihel,153), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1, -1/
  data (nhel(ihel,154), ihel = 1, 8) / 1, -1, -1,  1,  1, -1, -1,  1/
  data (nhel(ihel,155), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1, -1/
  data (nhel(ihel,156), ihel = 1, 8) / 1, -1, -1,  1,  1, -1,  1,  1/
  data (nhel(ihel,157), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1, -1/
  data (nhel(ihel,158), ihel = 1, 8) / 1, -1, -1,  1,  1,  1, -1,  1/
  data (nhel(ihel,159), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1, -1/
  data (nhel(ihel,160), ihel = 1, 8) / 1, -1, -1,  1,  1,  1,  1,  1/
  data (nhel(ihel,161), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1, -1/
  data (nhel(ihel,162), ihel = 1, 8) / 1, -1,  1, -1, -1, -1, -1,  1/
  data (nhel(ihel,163), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1, -1/
  data (nhel(ihel,164), ihel = 1, 8) / 1, -1,  1, -1, -1, -1,  1,  1/
  data (nhel(ihel,165), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1, -1/
  data (nhel(ihel,166), ihel = 1, 8) / 1, -1,  1, -1, -1,  1, -1,  1/
  data (nhel(ihel,167), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1, -1/
  data (nhel(ihel,168), ihel = 1, 8) / 1, -1,  1, -1, -1,  1,  1,  1/
  data (nhel(ihel,169), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1, -1/
  data (nhel(ihel,170), ihel = 1, 8) / 1, -1,  1, -1,  1, -1, -1,  1/
  data (nhel(ihel,171), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1, -1/
  data (nhel(ihel,172), ihel = 1, 8) / 1, -1,  1, -1,  1, -1,  1,  1/
  data (nhel(ihel,173), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1, -1/
  data (nhel(ihel,174), ihel = 1, 8) / 1, -1,  1, -1,  1,  1, -1,  1/
  data (nhel(ihel,175), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1, -1/
  data (nhel(ihel,176), ihel = 1, 8) / 1, -1,  1, -1,  1,  1,  1,  1/
  data (nhel(ihel,177), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1, -1/
  data (nhel(ihel,178), ihel = 1, 8) / 1, -1,  1,  1, -1, -1, -1,  1/
  data (nhel(ihel,179), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1, -1/
  data (nhel(ihel,180), ihel = 1, 8) / 1, -1,  1,  1, -1, -1,  1,  1/
  data (nhel(ihel,181), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1, -1/
  data (nhel(ihel,182), ihel = 1, 8) / 1, -1,  1,  1, -1,  1, -1,  1/
  data (nhel(ihel,183), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1, -1/
  data (nhel(ihel,184), ihel = 1, 8) / 1, -1,  1,  1, -1,  1,  1,  1/
  data (nhel(ihel,185), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1, -1/
  data (nhel(ihel,186), ihel = 1, 8) / 1, -1,  1,  1,  1, -1, -1,  1/
  data (nhel(ihel,187), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1, -1/
  data (nhel(ihel,188), ihel = 1, 8) / 1, -1,  1,  1,  1, -1,  1,  1/
  data (nhel(ihel,189), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1, -1/
  data (nhel(ihel,190), ihel = 1, 8) / 1, -1,  1,  1,  1,  1, -1,  1/
  data (nhel(ihel,191), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1, -1/
  data (nhel(ihel,192), ihel = 1, 8) / 1, -1,  1,  1,  1,  1,  1,  1/
  data (nhel(ihel,193), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1, -1/
  data (nhel(ihel,194), ihel = 1, 8) / 1,  1, -1, -1, -1, -1, -1,  1/
  data (nhel(ihel,195), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1, -1/
  data (nhel(ihel,196), ihel = 1, 8) / 1,  1, -1, -1, -1, -1,  1,  1/
  data (nhel(ihel,197), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1, -1/
  data (nhel(ihel,198), ihel = 1, 8) / 1,  1, -1, -1, -1,  1, -1,  1/
  data (nhel(ihel,199), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1, -1/
  data (nhel(ihel,200), ihel = 1, 8) / 1,  1, -1, -1, -1,  1,  1,  1/
  data (nhel(ihel,201), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1, -1/
  data (nhel(ihel,202), ihel = 1, 8) / 1,  1, -1, -1,  1, -1, -1,  1/
  data (nhel(ihel,203), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1, -1/
  data (nhel(ihel,204), ihel = 1, 8) / 1,  1, -1, -1,  1, -1,  1,  1/
  data (nhel(ihel,205), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1, -1/
  data (nhel(ihel,206), ihel = 1, 8) / 1,  1, -1, -1,  1,  1, -1,  1/
  data (nhel(ihel,207), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1, -1/
  data (nhel(ihel,208), ihel = 1, 8) / 1,  1, -1, -1,  1,  1,  1,  1/
  data (nhel(ihel,209), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1, -1/
  data (nhel(ihel,210), ihel = 1, 8) / 1,  1, -1,  1, -1, -1, -1,  1/
  data (nhel(ihel,211), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1, -1/
  data (nhel(ihel,212), ihel = 1, 8) / 1,  1, -1,  1, -1, -1,  1,  1/
  data (nhel(ihel,213), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1, -1/
  data (nhel(ihel,214), ihel = 1, 8) / 1,  1, -1,  1, -1,  1, -1,  1/
  data (nhel(ihel,215), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1, -1/
  data (nhel(ihel,216), ihel = 1, 8) / 1,  1, -1,  1, -1,  1,  1,  1/
  data (nhel(ihel,217), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1, -1/
  data (nhel(ihel,218), ihel = 1, 8) / 1,  1, -1,  1,  1, -1, -1,  1/
  data (nhel(ihel,219), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1, -1/
  data (nhel(ihel,220), ihel = 1, 8) / 1,  1, -1,  1,  1, -1,  1,  1/
  data (nhel(ihel,221), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1, -1/
  data (nhel(ihel,222), ihel = 1, 8) / 1,  1, -1,  1,  1,  1, -1,  1/
  data (nhel(ihel,223), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1, -1/
  data (nhel(ihel,224), ihel = 1, 8) / 1,  1, -1,  1,  1,  1,  1,  1/
  data (nhel(ihel,225), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1, -1/
  data (nhel(ihel,226), ihel = 1, 8) / 1,  1,  1, -1, -1, -1, -1,  1/
  data (nhel(ihel,227), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1, -1/
  data (nhel(ihel,228), ihel = 1, 8) / 1,  1,  1, -1, -1, -1,  1,  1/
  data (nhel(ihel,229), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1, -1/
  data (nhel(ihel,230), ihel = 1, 8) / 1,  1,  1, -1, -1,  1, -1,  1/
  data (nhel(ihel,231), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1, -1/
  data (nhel(ihel,232), ihel = 1, 8) / 1,  1,  1, -1, -1,  1,  1,  1/
  data (nhel(ihel,233), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1, -1/
  data (nhel(ihel,234), ihel = 1, 8) / 1,  1,  1, -1,  1, -1, -1,  1/
  data (nhel(ihel,235), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1, -1/
  data (nhel(ihel,236), ihel = 1, 8) / 1,  1,  1, -1,  1, -1,  1,  1/
  data (nhel(ihel,237), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1, -1/
  data (nhel(ihel,238), ihel = 1, 8) / 1,  1,  1, -1,  1,  1, -1,  1/
  data (nhel(ihel,239), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1, -1/
  data (nhel(ihel,240), ihel = 1, 8) / 1,  1,  1, -1,  1,  1,  1,  1/
  data (nhel(ihel,241), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1, -1/
  data (nhel(ihel,242), ihel = 1, 8) / 1,  1,  1,  1, -1, -1, -1,  1/
  data (nhel(ihel,243), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1, -1/
  data (nhel(ihel,244), ihel = 1, 8) / 1,  1,  1,  1, -1, -1,  1,  1/
  data (nhel(ihel,245), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1, -1/
  data (nhel(ihel,246), ihel = 1, 8) / 1,  1,  1,  1, -1,  1, -1,  1/
  data (nhel(ihel,247), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1, -1/
  data (nhel(ihel,248), ihel = 1, 8) / 1,  1,  1,  1, -1,  1,  1,  1/
  data (nhel(ihel,249), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1, -1/
  data (nhel(ihel,250), ihel = 1, 8) / 1,  1,  1,  1,  1, -1, -1,  1/
  data (nhel(ihel,251), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1, -1/
  data (nhel(ihel,252), ihel = 1, 8) / 1,  1,  1,  1,  1, -1,  1,  1/
  data (nhel(ihel,253), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1, -1/
  data (nhel(ihel,254), ihel = 1, 8) / 1,  1,  1,  1,  1,  1, -1,  1/
  data (nhel(ihel,255), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1, -1/
  data (nhel(ihel,256), ihel = 1, 8) / 1,  1,  1,  1,  1,  1,  1,  1/
  sqq_tt_bbeevv_ew = 0d0
  ntry = ntry + 1
  do ihel = 1, ncomb
    if (goodhel(ihel) .or. ntry < 10) then
      t = qq_tt_bbeevv_ew(iq, jf, p1, p2, p3, p4, p5, p6, p7, p8, nhel(1, ihel), channel)
      sqq_tt_bbeevv_ew = sqq_tt_bbeevv_ew + t
      if (t > 0d0 .and. .not. goodhel(ihel)) then
          goodhel(ihel)= .true.
      endif
    end if
  enddo
  sqq_tt_bbeevv_ew = sqq_tt_bbeevv_ew / 4.d0
end function sqq_tt_bbeevv_ew


function qq_tt_bbeevv_ew(iq, jf, p1, p2, p3, p4, p5, p6, p7, p8, nhel, channel)

  ! function generated by madgraph
  ! returns amplitude squared summed/avg over colors
  ! for the point in phase space p1, p2, p3, p4, p5, p6, p7, p8
  ! and helicity nhel(1), nhel(2)
  ! for process: q q~ -> A, Z, Z' -> b b~ e+ e- v v~

  use kinds
  use configuration, only: include_a, include_z, include_x, interference
  use modelling

  implicit none

  ! constants
  integer, parameter :: ngraphs = 7, nexternal = 8
  real, parameter :: zero = 0d0

  ! arguments
  integer :: iq, jf
  real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
  integer :: nhel(nexternal)
  integer, intent(in), optional :: channel

  ! variables
  integer :: i, j
  complex(kind=complex) amp_tmp, amp_tmp2
  complex(kind=complex) amp(ngraphs)
  complex(kind=complex) w1(6), w2(6), w3(6), w4(6), w5(6)
  complex(kind=complex) w6(6), w7(6), w8(6), w9(6), w10(6)
  complex(kind=complex) w11(6), w12(6), w13(6), w14(6), w15(6)
  real(kind=default) :: gAq(2)
  real(kind=default) :: gZq(2)
  real(kind=default) :: gxq(2,5)
  real(kind=default) :: gxq_tmp(2), gxf_tmp(2)
  real(kind=default) :: atot2, asum2

  ! initial state

  ! uu or cc
  if ((iq == 3) .or. (iq == 7)) then
    do i = 1, 2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1, 5
        gxq(i,j) = gxu(i,j)
      end do
    enddo

  ! dd or ss
  else if ((iq == 4) .or. (iq == 8)) then
    do i = 1, 2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1, 5
        gxq(i,j) = gxd(i,j)
      end do
    enddo

  ! tt
  else if (iq == 11) then
    do i = 1, 2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1, 5
        gxq(i,j) = gxt(i,j)
      end do
    enddo

  ! bb
  else if (iq == 12) then
    do i = 1, 2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1, 5
        gxq(i,j) = gxb(i,j)
      end do
    enddo
  end if

  ! initialise amplitudes
  do i = 1, ngraphs
    amp(i) = 0.d0
  enddo

  ! wavefunctions
  call ixxxxx(p1, fmass(iq), nhel(1),  1, w1)
  call oxxxxx(p2, fmass(iq), nhel(2), -1, w2)
  call oxxxxx(p3, fmass(12), nhel(3),  1, w3)
  call ixxxxx(p4, fmass(12), nhel(4), -1, w4)
  call ixxxxx(p5, fmass(1), nhel(5), -1, w5)
  call oxxxxx(p6, fmass(2), nhel(6),  1, w6)
  call oxxxxx(p7, fmass(1), nhel(7),  1, w7)
  call ixxxxx(p8, fmass(2), nhel(8), -1, w8)

  ! W coupled to e- nu~ vector current
  call jioxxx(w5, w7, gWf, wmass, wwidth, w10)

  ! W coupled to e+ nu vector current
  call jioxxx(w8, w6, gWf, wmass, wwidth, w11)

  ! t coupled to b W vector current
  call fvoxxx(w3, w10, gWf, fmass(11), fwidth(11), w12)

  ! t~ coupled to b~ and W vector current
  call fvixxx(w4, w11, gWf, fmass(11), fwidth(11), w13)

  if (include_a) then
    ! A coupled to qq~ vector current
    call jioxxx(w1, w2, gAq, amass, amass, w9)
    ! A diagram
    call iovxxx(w13, w12, w9, gAu, amp(1))
  end if

  if (include_z) then
    ! Z coupled to qq~ vector current
    call jioxxx(w1, w2, gZq, zmass, zwidth, w14)
    ! Z diagram
    call iovxxx(w13, w12, w14, gZu, amp(2))
  end if

  ! Z' diagrams
  if (include_x) then
    do i = 1, 5
      if (xmass(i) > 0) then
        do j = 1, 2
          gxq_tmp(j) = gxq(j,i)
          gxf_tmp(j) = gxt(j,i)
        end do
        call jioxxx(w1, w2, gxq_tmp, xmass(i), xwidth(i), w15)
        call iovxxx(w13, w12, w15, gxf_tmp, amp(2 + i))
      end if
    end do
  end if

  ! total M*M for given helicity combination
  atot2 = 0.d0
  amp_tmp = (0.d0, 0.d0)
  amp_tmp2 = (0.d0, 0.d0)

  if (interference == 0) then
    do i = 1, ngraphs
      atot2 = atot2 + amp(i) * conjg(amp(i))
    end do
    qq_tt_bbeevv_ew = atot2

  else if (interference == 1) then
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    atot2 = atot2 + amp_tmp * conjg(amp_tmp)

    if (present(channel) .and. channel > 0) then
        asum2 = 0.d0
        do i = 1, ngraphs
          asum2 = asum2 + amp(i) * conjg(amp(i))
        end do

        if (asum2 > 0) then
            qq_tt_bbeevv_ew = amp(channel) * conjg(amp(channel)) * atot2 / asum2
        else
            qq_tt_bbeevv_ew = 0
        end if
    else
        qq_tt_bbeevv_ew = atot2
    end if

    else if (interference == 2) then
        do i = 1, 2
          amp_tmp = amp_tmp + amp(i)
        end do
        atot2 = atot2 + amp_tmp * conjg(amp_tmp)
        do i = 3, ngraphs
          atot2 = atot2 + amp(i) * conjg(amp(i))
        end do
        qq_tt_bbeevv_ew = atot2

  else if (interference == 3) then
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    atot2 = atot2 + amp_tmp * conjg(amp_tmp) - amp_tmp2 * conjg(amp_tmp2)
    qq_tt_bbeevv_ew = atot2

  else if (interference == 4) then
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    atot2 = atot2 + amp_tmp * conjg(amp_tmp) - amp_tmp2 * conjg(amp_tmp2)
    do i = 3, ngraphs
      atot2 = atot2 - amp(i) * conjg(amp(i))
    end do
    qq_tt_bbeevv_ew = atot2
  end if


end function qq_tt_bbeevv_ew

end module tt_bbeevv
