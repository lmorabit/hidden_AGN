rule testplot:
	input:
		"src/static/RLFS_50MYR_SF.csv"
		"src/static/RLFS_50MYR_AGN.csv"
	output:
		"test.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/testplot.py"
rule mauchsadler:
	input:
		"src/static/mauch_sadler_table5.csv"
		"src/static/lockman_final_cross_match_catalogue-v1.0.fits"
		"src/static/lockman_rms_starmask_optical.fits"
	output:
		"mauch_sadler_RLFs.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/mauchsadler.py"
