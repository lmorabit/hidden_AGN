rule generateZmaxes:
	input:
		"src/static/lockman_final_cross_match_catalogue-v1.0.fits"
		"src/static/lockman_rms_starmask_optical.fits"
	output:
		"src/data/zmaxes.fits"
	conda:
		"environment.yml"
	script:
		"src/scripts/generateZmaxes.py"
