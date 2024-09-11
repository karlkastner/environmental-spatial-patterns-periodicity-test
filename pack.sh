#!/bin/bash
d=$(basename $(pwd))
cd ../
zip  -9 $d/$d.zip -r $d \
	--exclude \
		$d'/.git/*' \
		$d'/.gitignore' \
		$d'/img/*' \
		$d'/other/*' \
		$d'/output/*.mat' \
                $d'/mat/*'

