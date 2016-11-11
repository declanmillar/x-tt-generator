function sgg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8)

  ! returns amplitude squared summed/avg over colors and helicities
  ! for the point in phase space p1, p2, p3, p4, ...
  ! for process: g g -> t t~ -> b b~ ta+ ta- vt vt~

  implicit none

  ! functions
  real*8 :: sgg_tt_bbeevv
  real*8 :: gg_tt_bbeevv

  ! constants
  integer, parameter :: nexternal = 8, ncomb = 256

  ! arguments
  real*8 :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

  ! local variables
  integer :: i, j
  integer :: nhel(nexternal, ncomb), ntry
  real*8 :: t
  integer :: ihel
  logical :: goodhel(ncomb)
  data goodhel /ncomb*.false./
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
      t = gg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8, nhel(1, ihel))
      sgg_tt_bbeevv = sgg_tt_bbeevv + t
      if (t > 0.d0 .and. .not. goodhel(ihel)) then
          goodhel(ihel)= .true.
      endif
    end if
  enddo
  sgg_tt_bbeevv = sgg_tt_bbeevv / 4.d0
  ! if (sgg_tt_bbeevv == 0.d0) return
end function sgg_tt_bbeevv


function gg_tt_bbeevv(p1, p2, p3, p4, p5, p6, p7, p8, nhel)

  ! returns amplitude squared summed/avg over colors
  ! for the point in phase space p1, p2, p3, p4, ...
  ! and helicity nhel(1), nhel(2), ...
  ! for process : g g  -> t t~ -> b b~ ta+ ta- vt vt~

  use modelling

  implicit none

  real*8 :: gg_tt_bbeevv

  ! constants
  integer, parameter :: ngraphs = 3, neigen = 2, nexternal = 8
  real*8, parameter :: zero = 0.d0

  ! arguments
  real*8 :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)
  integer :: nhel(nexternal)

  ! local variables
  integer :: i, j
  real*8 :: eigen_val(neigen), eigen_vec(ngraphs, neigen)
  complex*16 ztemp
  complex*16 amp(ngraphs)
  complex*16 w1(6), w2(6), w3(6), w4(6), w5(6)
  complex*16 w6(6), w7(6), w8(6), w9(6), w10(6)
  complex*16 w11(6), w12(6), w13(6), w14(6), w15(6)

  ! initialise amplitudes
  do i = 1, ngraphs
    amp(i) = 0.d0
  enddo

  ! color data
  data eigen_val(1) /7.2916666666666588e-02/
  data eigen_vec(1, 1) /7.0710678118654768e-01/
  data eigen_vec(2, 1) /7.0710678118654735e-01/
  data eigen_vec(3, 1) /0.0000000000000000e+00/
  data eigen_val(2) /2.8125000000000000e-01/
  data eigen_vec(1, 2) /-4.0824829046386313e-01/
  data eigen_vec(2, 2) /4.0824829046386285e-01/
  data eigen_vec(3, 2) /8.1649658092772603e-01/

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

  gg_tt_bbeevv = 0.d0
  do i = 1, neigen
    ztemp = (0.d0, 0.d0)
    do j = 1, ngraphs
      ztemp = ztemp + eigen_vec(j,i)*amp(j)
    enddo
    gg_tt_bbeevv = gg_tt_bbeevv + ztemp*eigen_val(i)*conjg(ztemp)
  enddo
  ! call gaugecheck(amp, ztemp, eigen_vec, eigen_val, ngraphs, neigen)
end function gg_tt_bbeevv
