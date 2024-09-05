#!/bin/bash
d=$(basename $(pwd))
cd ../
zip  -9 $d/$d.zip -r $d \
	--exclude $d'/mat/*' $d'/other/*' $d'/.git/*' $d'/img/*'

