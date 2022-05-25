**Zadatak**. Zadana je domena ![Omega](./doc/o.png)
 (datoteka `transfinite-3D-simplex.msh`) prikazana na slici:

![domna](./doc/domena.png)


Riješavamo zadaću distribucije temperature  domeni ![Omega](./doc/o.png).
Materijal je homogen i nema vanjskog volumnog izvora topine,  tako da
temperatura `u` zadovoljava Laplaceouvu jednadžbu. Na gornjoj plohi
granice ![pOmega](./doc/do.png), `z=1`, zadana je stalna  temperatura od 100 stupnjeva C
(temperaturna skala je Celzijeva). Bočne stranice su termički izolirane,
a donja zakrivljena strana zrači toplinu prema Stefan-Boltzmanovom 
zakonu, tako da na njoj imamo rubni uvjet:

![rubni uvjet](./doc/eqn.png).

Konstanta c iznosi 1E-8, a ambijentalna temperatura ![ue](./doc/ue.png) iznosi 25 stupnjeva C.
Izračunajte distribuciju temperature u domeni  ![Omega](./doc/o.png).


