;  25/April/2013   J.A.Rubino
;
;


FUNCTION pars2pmodel_structure, a
  pmodel = { P_f:a[0], P_x:a[1], P_y:a[2], P_c:a[3], P_n:a[4], P_a:a[5], P_b:a[6] } 
  return,pmodel
END

;-----------------

PRO sancho_mfi_pmodel, pmodel, ihorn=ihorn, center=center, zero=zero

; Keywords
if not keyword_set(ihorn) then ihorn=0
if not keyword_set(center) then center=0
if not keyword_set(zero) then zero=0


; Initializing the structure to zero.
get_params_pointing_model, pmodel

arsec2rad  = !DPI/3600.0/180.0 

; Datos de Marta (email 17/Ene/2013)
;  IA    +3606.98
;  IE      +34.63
;  NPAE    -45.16
;  CA    -2942.01
;  AN      +39.62
;  AW      -81.51
;  TF      -38.45
pmodel.P_f =   -38.45 * arsec2rad
pmodel.P_x =    39.62 * arsec2rad
pmodel.P_y =   -81.51 * arsec2rad
pmodel.P_c = -2942.01 * arsec2rad
pmodel.P_n =   -45.16 * arsec2rad
pmodel.P_a =  3606.98 * arsec2rad
pmodel.P_b =    34.63 * arsec2rad


; Solucion de Denis (25/April/2013)
arsec2rad  = !DPI/3600.0/180.0 
pmodel.P_f = -2732.71 * arsec2rad
pmodel.P_x =   571.74 * arsec2rad
pmodel.P_y =   546.24 * arsec2rad
pmodel.P_c =  -621.24 * arsec2rad
pmodel.P_n = -1355.04 * arsec2rad
pmodel.P_a =   573.32 * arsec2rad
pmodel.P_b = -2521.11 * arsec2rad


; Solucion de prueba aproximada (25/Abril/2013)
arsec2rad  = !DPI/3600.0/180.0 
pmodel.P_f =    0.0 * arsec2rad
pmodel.P_x =    0.0 * arsec2rad
pmodel.P_y =    0.0 * arsec2rad
pmodel.P_c =    0.0 * arsec2rad
pmodel.P_n =    0.0  * arsec2rad
pmodel.P_a =    0.0 * arsec2rad
pmodel.P_b = -1500.0 * arsec2rad


; Solucion de Denis (16/Dec/2013), ajustando CASS y CRAB
; para la bocina1.
pmodel.P_f =   498.119049 * arsec2rad
pmodel.P_x =  -155.344711 * arsec2rad
pmodel.P_y =     3.690629 * arsec2rad
pmodel.P_c =  2086.427246 * arsec2rad
pmodel.P_n = -2668.219727 * arsec2rad
pmodel.P_a = -1046.779175 * arsec2rad
pmodel.P_b = -1191.227661 * arsec2rad


; Solucion individual por bocinas de Denis. 17/Dec/2013. Habia error de JD.
horn1 = [505.197937, -157.148483,   0.000000, 2064.570312, -2651.906006, -1029.397095, -1188.395142 ] * arsec2rad
horn2 = [ 72.298409,  -37.873341, -62.718033, 1911.634277, -2102.814697,  -678.208862, -1733.438965 ] * arsec2rad
horn3 = [153.467957, -163.480942, -83.365723,   68.297401,  -703.291077,   652.511536, -1453.635864 ] * arsec2rad
horn4 = [464.715668, -167.701447, -34.419003,    0.000000,  -819.801331,   520.852905, -1053.586548 ] * arsec2rad


; Nueva solucion por bocinas, combinando CRAB y CASS. Denis, 28/Dec/2013. Sin
; error de JD
horn1 = [341.537476, -110.144310, 104.277802, 484.378876, -1253.372559,   0.000000, -1270.996338 ] * arsec2rad
horn2 = [ 53.119659,  -54.753948,   0.000000,   0.000000,  -458.266998, 591.758179, -1764.614258 ] * arsec2rad
horn3 = [ 30.475422,  -42.224468, -15.681525,   0.000000,  -325.047363, 508.562042, -1485.779541 ] * arsec2rad
horn4 = [289.966553,  -30.458158,   7.766873, 492.643372,  -885.281860,   0.000000, -1121.425415 ] * arsec2rad

; Nueva solucion por bocinas, usando solo CRAB. Denis, 28/Dec/2013. Sin error de
; JD 
horn1 = [464.396729, -370.451538,  59.979919, 1782.341431, -2477.275879, -1068.910645, -1081.417969 ] * arsec2rad
horn2 = [335.435333, -209.219330, -76.672783,  739.692749, -1166.342163,     0.000000, -1633.830200 ] * arsec2rad
horn3 = [188.409836,  -78.538559, -76.420784, 1049.668213, -1287.030396,  -160.660873, -1403.562988 ] * arsec2rad
horn4 = [495.030884, -191.893875, -34.547272,  660.684204, -1228.798828,     0.000000,  -971.630310 ] * arsec2rad


; Calculo por pasos: P_b --> P_a --> P_f --> all. Usando CRAB+CASS. Denis, 8/Ene/2014
horn1 = [ 293.954478, -77.971687, 120.967911, 977.698057, -1605.863941, -379.206695, -1290.334907 ] *  arsec2rad
horn2 = [  -9.048263, -39.855222,  18.266917, 506.363741,  -822.873307,  199.387915, -1798.069416 ] *  arsec2rad
horn3 = [  15.643603, -29.163693,  -5.156503, 356.247605,  -582.044545,  234.751479, -1488.186947 ] *  arsec2rad
horn4 = [ 267.871644, -12.544701,  14.878328, 598.688323,  -961.246133,  -78.915599, -1126.947130 ] *  arsec2rad


; Nuevo modelo. Denis, 13/Ene/2014. Corregido error precision en JD. El mejor
; hasta ahora. Usa focal=3637mm. 
horn1 = [ 361.197313, 81.906147, -27.720241, -3891.063433, 891.804788,  4509.580503, -1225.131212 ] *  arsec2rad
horn2 = [-271.864054, 56.644465, -32.972026, -1935.401913,   0.000000,  3108.027097, -1978.780716 ] *  arsec2rad        
horn3 = [-157.777922, 62.584273, -39.353865, -2128.438709,   0.000000,  3347.740267, -1590.406421 ] *  arsec2rad        
horn4 = [ 534.826868, 64.082589, -29.267457, -2494.121632,   0.000000,  3602.686003,  -915.940276 ] *  arsec2rad 


; Nuevo modelo. Fit simultaneo a las 4 bocinas. Denis. 15/Enero/2014. Malo para
; algunas bocinas.
horn1 = [ 242.571640, 80.881378, -31.033945, -2726.602783,   0.000000,  3927.132812, -1302.101807 ] *  arsec2rad 
horn2 = horn1
horn3 = horn1
horn4 = horn1


; Fit siguiendo el esquema de TPOINT. 23/Enero/2014. Malo. Offset y elipticidad
horn1 = [ 24.817992, -61.836121, 91.388830, -94.221689,	-249.332040, -69.433542, -1366.572226 ] *  arsec2rad
horn2 = [ 24.818005, -61.836155, 91.388879, 107.071513,	-249.332173, -69.433579, -1737.406523 ] *  arsec2rad
horn3 = [ 24.818003, -61.836147, 91.388867, 232.602901,	-249.332142, -69.433571, -1666.289745 ] *  arsec2rad
horn4 = [ 24.817995, -61.836130, 91.388842, 221.753339,	-249.332070, -69.433551, -1379.688600 ] *  arsec2rad


; Fit por bocinas con nueva focal de 3700mm. 27/Ene/2014.
horn1 = [ 319.146146, 83.151900, -23.616893, -2716.493422,   0.000000, 3852.068745, -1112.758258 ] *  arsec2rad
horn2 = [-286.716442, 56.604657, -34.769469, -2103.255488,   0.000000, 3125.197762, -1908.697116 ] *  arsec2rad
horn3 = [   0.000000, 52.629357, -50.119368, -3219.108120, 776.883238, 4017.435764, -1634.993920 ] *  arsec2rad
horn4 = [ 514.127830, 62.949540, -30.353121, -2358.211731,   0.000000, 3618.475458, -1027.294733 ] *  arsec2rad



; Fit usando esquema de TPOINT, dejando Pc, Pn y Pb depender de la bocina. Focal
; 3700mm. 27/Enero/2014
horn1 = [ 24.819039, -61.838731, 91.392685, 3344.882664, -4816.837694, -69.436471, -1359.189363 ] *  arsec2rad
horn2 = [ 24.819296, -61.839370, 91.393630, 2604.863526, -3610.772011, -69.437189, -1732.308296 ] *  arsec2rad
horn3 = [ 24.819302, -61.839385, 91.393653, 2498.378652, -3394.819605, -69.437206, -1659.790492 ] *  arsec2rad
horn4 = [ 24.818780, -61.838085, 91.391731, 2878.540468, -3889.527842, -69.435746, -1379.467753 ] *  arsec2rad


; Fit 5/Feb/2014, usando esquema TPOINT, dejando Pc, Pn, Pb libres. Focal
; 3700mm. Denis.
horn1 = [ 25.017630, -26.427459, 46.129333, 3291.258680, -4552.873449, -174.404276, -1343.539231 ] *  arsec2rad
horn2 = [ 25.017986, -26.427835, 46.129989, 3378.249552, -4527.936296, -174.406758, -1788.061742 ] *  arsec2rad
horn3 = [ 25.018917, -26.428818, 46.131706, 3515.977221, -4691.506930, -174.413246, -1719.600660 ] *  arsec2rad
horn4 = [ 25.019724, -26.429670, 46.133194, 3667.568499, -4878.251937, -174.418872, -1377.741084 ] *  arsec2rad


; Fit 5/Feb/2014 individual por bocinas. Todos los datos (incluidos
; nominal). Denis. Focal 3700mm. 
; This is currently the best solution.
horn1 = [ 645.055552, 110.473721, -78.897113, -2437.223362, 0.000000, 3519.963935,  -891.423134 ] *  arsec2rad
horn2 = [-559.842110,  56.687210, -49.288305, 0.000000, -1721.153046, 1870.611004, -2108.336180 ] *  arsec2rad 
horn3 = [-379.809842,  46.793957, -61.435784, -2310.235444, 0.000000, 3521.138938, -1925.361772 ] *  arsec2rad
horn4 = [ 565.595688,  73.590737, -20.853737, -2364.792847, 0.000000, 3620.518854,  -987.025009 ] *  arsec2rad


; Fit 7/dic/2016. Solucion por bocinas, pero transformando antes el centro
; del plano focal. Mejor solucion nueva, desde Dic/2016
if (center) then begin
	horn1 = [  719.26, 115.44, -94.45,  -241.70, -2441.73, 2649.56,  -894.51 ] *  arsec2rad
	horn2 = [ -524.86,  52.59, -51.57, -2958.06,   104.49, 4068.00, -2142.27 ] *  arsec2rad
	horn3 = [ -483.52,  36.70, -87.57, -1907.00,    59.82, 2897.05, -1960.48 ] *  arsec2rad
	horn4 = [  728.59,  49.57, -82.04, -2797.36,   631.80, 3740.33,  -938.26 ] *  arsec2rad
endif


;-------------------------------------------
; Arrange pmodel in the correct format
pmodel1 = pars2pmodel_structure(horn1)
pmodel2 = pars2pmodel_structure(horn2)
pmodel3 = pars2pmodel_structure(horn3)
pmodel4 = pars2pmodel_structure(horn4)

txt = "     > Using pointing model for horn number "+string(ihorn,format='(i1)')
case ihorn of 
   1: pmodel = pmodel1 
   2: pmodel = pmodel2 
   3: pmodel = pmodel3 
   4: pmodel = pmodel4 
   else: txt = "     > Using default pointing model"
endcase


;==========================
; Set pmodel to zero if requested
if zero then begin
	get_params_pointing_model, pmodel
	txt =  "     > Using null pointing model (all parameters to zero)."
endif
print,txt


END
