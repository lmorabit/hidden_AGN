rule testplot:
	input:
		"RLFS_50MYR_SF.csv"
		"RLFS_50MYR_AGN.csv"
		"mauch_sadler_table5.csv"
	output:
		"test.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/testplot.py"
