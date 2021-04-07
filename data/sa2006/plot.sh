!#/bin/ bash
cd ./gmts/
gmt begin plot pdf
  gmt coast -Rg -Ggrey -Sazure -B -JR15c -Wthin
  gmt psxy alaska.aleutians.contours.gmt -JR15c -Wblue
  gmt psxy andes.contours.gmt -JR15c -Wblue
  gmt psxy central.america.contours.gmt -JR15c -Wblue
  gmt psxy kamchatka.marianas.contours.gmt -JR15c -Wblue
  gmt psxy kyushu.ryukyu.contours.gmt -JR15c -Wblue
  gmt psxy lesser.antilles.contours.gmt -JR15c -Wblue
  gmt psxy n.philippines.contours.gmt -JR15c -Wblue
  gmt psxy new.britain.solomon.contours.gmt -JR15c -Wblue
  gmt psxy s.philippines.contours.gmt -JR15c -Wblue
  gmt psxy scotia.contours.gmt -JR15c -Wblue
  gmt psxy sumatra.banda.sea.contours.gmt -JR15c -Wblue
  gmt psxy tonga.new.zealand.contours.gmt -JR15c -Wblue
  gmt psxy vanuatu.contours.gmt -JR15c -Wblue
gmt end show
