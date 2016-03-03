* Encoding: UTF-8.
* 0 41 30 79 67 87 80 149 111 212 234 213 232 258 286 274.
* item 10, item 3, item 8, item 1, item 16, item 13, item 15, item 7, item 14, item 2, item 11, item 12, item 5, item 4, item 9, item 6.

compute controlT=
  V7*sin(171.5*r) + V14*sin(133.5*r) + V1*sin(101.5*r) + V15*sin(102.5*r) + V13*sin(109.5*r) + V16*sin(99.5*r) +
  V3*sin(63.5*r) + V8*sin(52.5*r) + V6*sin(296.5*r) + V10*sin(22.5*r) + V4*sin(280.5*r) + V9*sin(308.5*r) +
  V5*sin(254.5*r) + V11*sin(256.5*r) + V2*sin(234.5*r) + V12*sin(235.5*r).
VARIABLE LABELS controlT 'control score leerkrachten gebasseerd op theoretische hoekwaardes'.
execute.

compute affiliationT=
 V7*cos(171.5*r) + V14*cos(133.5*r) + V1*cos(101.5*r) + V15*cos(102.5*r) + V13*cos(109.5*r) + V16*cos(99.5*r) +
  V3*cos(63.5*r) + V8*cos(52.5*r) + V6*cos(296.5*r) + V10*cos(22.5*r) + V4*cos(280.5*r) + V9*cos(308.5*r) +
  V5*cos(254.5*r) + V11*cos(256.5*r) + V2*cos(234.5*r) + V12*cos(235.5*r).
Variable labels affiliationT 'affiliation score per leerkracht gebaseerd op theoretische hoekwaarden'.
execute. 

* hoeken tussen 0 en 90 graden (DUS: y = positief en x = positief). omrekenen van radialen naar graden = * 57.298.
do if controlT GT 0 AND affiliationT GT 0.
compute Person.score = 57.298 * ARTAN(controlT / affiliationT).
end if. 

do if controlT LT 0 AND affiliationT GT 0.
compute Person.score = 57.298 * ((2*3.1416) + ARTAN(controlT / affiliationT)).
end if. 

do if controlT GT 0 AND affiliationT LT 0.
compute Person.score = 57.298 * (3.1416 + ARTAN(controlT / affiliationT)).
end if. 

do if controlT LT 0 AND affiliationT LT 0.
compute Person.score = 57.298 * (3.1416 + ARTAN(controlT / affiliationT)).
end if. 
