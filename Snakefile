rule testplot:
	input:
		"src/static/RLFS_50MYR_SF.csv"
		"src/static/RLFS_50MYR_AGN.csv"
		"src/static/mauch_sadler_table5.csv"
	output:
		"test.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/testplot.py"
