#!/usr/bin/perl

my $start;
$start = time();

print "GA - START\n";

# Parameters
my $iterationnumber=25;		# number of iterations
my $crossoverprob=0.3;		# crossover probability
my $mutationprob=0.025;		# mutation probability
my $fitnesspow=3.0;		# pow value for fitness formula
my $indnumber=20;		# number of individuals (solutions)
my $holdoutval=0.33333;		# percentage of inputs that will be used for cross validation
my $infile="./input.txt";	# indicators results input file

# Variables
my $linecount=-1;
my $maxval=0;
my $end;
my $indsize;			# Size of individuals - number of chromossomes
my $crossresultssize=0;
my @ind;
my @indaux;
my @results;
my @fields;
my @crossresults;
my @crossval;
my @loadedresults;
my @results1;

use Time::HiRes qw( time );
use POSIX qw/ceil/;

#system(clear);

# Reading parameters
# Reading input file
if ( $ARGV[0] ){
	$infile = $ARGV[0];
}
# Reading population size
if ( $ARGV[1] ){
	$indnumber = $ARGV[1];
}

if ( $ARGV[2] ){
	$iterationnumber = $ARGV[2];
}
if ( $ARGV[3] ){
	$thingyflag = 1;
}

# Printing GA conditions
sub printgaconditions {
	print ">>> GA conditions\n";
	print "Number of iterations:\t\t$iterationnumber\n";
	print "Crossover probability:\t\t$crossoverprob\n";
	print "Mutation probability:\t\t$mutationprob\n";
	print "Holdout probability:\t\t$holdoutval\n";
	print "Fitness pow:\t\t\t$fitnesspow\n";
	print "Population size:\t\t$indnumber\n";
	print "Number of chromossomes:\t\t$indsize\n";
	print "Input file:\t\t\t$infile\n";
	$resultssize = scalar(@results);
	print "Training set size:\t\t$resultssize\n";
	$crossresultssize = scalar(@crossresults);
	print "Cross validation set size:\t$crossresultssize\n";
	print "\n";
}

# populates a matrix with the indicator results calculated in the spreadsheet
sub popresults {

	open(INFILE,$infile) or die("Cant open file:$!");
	while (<INFILE>){

	        $line=$_;
		chomp($line);
		@fields = split ',', $line;
#		print "fields: @fields\n";

		$linecount=$linecount+1;
#		print "linecount: $linecount\n";
	        for my $colcount (0..scalar(@fields)-1){
#			print "colcount: $colcount\n";
			$loadedresults[$linecount][$colcount] = @fields[$colcount];
		}
	}
	close(INFILE);

	# Defining individual size based on input file
	$indsize = scalar(@fields)-1;
	#print "indsize= $indsize\n";

}

# Holdout method to divide the loadedresults matrix in two parts
sub holdout {
	my $myrand=0.0;
	my $resultsrow=0;
	my $crossvalrow=0;
	for my $row (0..scalar(@loadedresults)-1){
		$myrand=rand();
		if ( $myrand gt $holdoutval ){
#			print "resultsrow :$resultsrow\n";
			$results[$resultsrow]=$loadedresults[$row];
			$resultsrow=$resultsrow+1;
		}
		else {
#			print "crossvalrow :$crossvalrow\n";
			$crossresults[$crossvalrow]=$loadedresults[$row];
			$crossvalrow=$crossvalrow+1;
		}
	}
	$crossresultssize = scalar(@crossresults)-1;
}

# Print arrays
sub printarray {

	my @array=@_;
	my $sum;

	for my $indn (0..scalar(@array)-1){
		$sum=0.0;
#		print "Line $indn: ";
		print "Line\t$indn:\t";
		for my $column (0..scalar(@{$array[0]})-1){
#			print "$array[$indn][$column] ";
			if ( $column >= $indsize ) {
#				print "| ";
				print "|\t";
			}
			printf "%.2f\t",$array[$indn][$column] ;
#			printf "%.2f ",$array[$indn][$column] ;
		}
#		print "| Sum: $sum\n";
		print "\n";
	}
	print "\n";
}

# Initialize subjects
sub initialize {

	my $value;
	for my $indn (0..$indnumber-1){
		for my $column (0..$indsize-1){

			my $myrand=rand();
			if ( $myrand gt 0.5 ){
				$value=1;
			}
			else {
				$value=0;
			}
			$ind[$indn][$column]=$value;
		}
		$ind[$indn][$indsize]=0.00;
		$ind[$indn][$indsize+1]=0.00;
	}
}

# Normalize fitness column for next generation selection through roulette wheel method
sub selectionprep {

	my $max=0.00;
	my $diffsum=0.00;
	my $fitnessval=0.00;
	my $diffabs=0.00;

	# Getting the max of fitness column
	for my $row (0..scalar(@ind)-1){
		$fitnessval=$ind[$row][$indsize];
		if ( $max < $fitnessval ) {
			$max=$fitnessval;
		}
	}
	# abs of diff between fitness and its max
	for my $row (0..scalar(@ind)-1){
		$diffabs=abs($ind[$row][$indsize]-$max);
		$ind[$row][($indsize+1)]=$diffabs ** $fitnesspow;;
		# Getting the sum of (fitness - max)
		$diffsum=$diffsum+$ind[$row][($indsize+1)];
	}

	return $diffsum;
}

# Fitness function
sub calcfitness {

	my $sum;
	my $sum2;
	my $max;

	# For each subject...
	for my $row (0..scalar(@ind)-1){
		$sum2=0.00;
		# For each result line...
		for my $resultline (0..(scalar(@results)-1)){
#			print "result line: $resultline\n";
			$sum=0.00;
			# Summation in n of r_m - (x_kn * y_mn) ### this is the fitness function
			for my $col (0..($indsize-1)){
				$sum=$sum + ($results[$resultline][scalar(@{$results[0]}) - 1] - ( $results[$resultline][$col] * $ind[$row][$col]));
			}
			# Summation in m of abs of variable sum (see second line above)
#			print "partial sum for individual $row, result line $resultline: $sum\n";
			$sum2=$sum2 + abs($sum);
		}
#		print "fitness for subject $row: $sum2\n";
		$ind[$row][$indsize]=$sum2;
#		print "total sum for ind $row: $sum2\n";
	}
#	print "\n";

	$max=&selectionprep;

	print ">>> Sorting array\n";
	&sortarray();

	return $max;
}


#Sort ind array by its last field (adjusted fitness)
sub sortarray {

	my @array = sort cmpfunc @ind;
	@ind=@array;

}
sub cmpfunc {
	my $size;
	$size=scalar(@{$ind[0]})-1;
	return (($a->[$size] <=> $b->[$size]) or 
		($a->[$size-1] <=> $b->[$size-1]));
}

# Select new population
sub selection {

	my $max=$_[0];
#	print "diffsum: $max\n";
	my $sum;

	# loop scalar(population) times to get the same number of subjects for next generation
	for my $count (0..scalar(@ind)-1){
		my $myrand=rand() * $max;
		my $myrand=int($myrand + 0.5);
#		print "# $count - rand: $myrand\n";

		$sum=0;
		# for each subject...
		for my $row (0..scalar(@ind)-1){
#			print "internal: $sum -- " . ($ind[$row][($indsize+1)] + $sum) . "\n";
			if ( ($myrand >= $sum ) and ($myrand <= ($ind[$row][($indsize+1)] + $sum) ) ) {
#				print "Choice $count: # $row\n";
				for my $col (0..($indsize+1)){
					$indaux[$count][$col]=$ind[$row][$col];
				}
				# If found, skip the rest of the loop
				last;
			}
			$sum=$sum+$ind[$row][($indsize+1)];
		}
	}
	@ind=@indaux;
	print "\n";
	print ">>> Sorting array\n";
	&sortarray();
}

# Crossover
sub crossover {
	my $pairs=0;
	my $aux;
	my $myrand=0;
	my $addone=0;
	$pairs=int($indnumber/2);
	# In case of odd number of individuals, add 50% chances of crosssing over pairs shifted by 1
	# ie., in case of 5 individuals, add 50% chances of crossing over paiasr (2-3) and (4-5)
	# without the code below, only pairs (1-2) and (3-4) would have a chance of crossing over
	$remainder=$indnumber%2;
#	print "remainder: $remainder\n";
	if ( $remainder ){
		$myrand=rand();
		if ( $myrand > 0.5 ){
			$addone=1;
		}
#		print "addone: $addone\n";
	}
#	print "crossoverprob: $crossoverprob\n";
#	print "# of pairs: $pairs\n";
	# For possible pairs loop...
	for my $count (1..$pairs) {
#		print "count: $count\n";

		# will crossover happen?
		$myrand=rand();
#		print "rand: $myrand\n";
		if ( $myrand <= $crossoverprob ) {

			$first=(2*$count)-1+$addone;
			$second=(2*$count)-2+$addone;
			my $myrandpos=ceil(rand()*($indsize-1));
#			my $myrandpos=ceil(rand()*4);
			print "Crossover between $first and $second at position $myrandpos\n";
#			print "myrandpos: $myrandpos\n";

			# From position pos till the last data position (scalar-3) swap values
			for my $pos ($myrandpos..($indsize-1)) {
				$aux=$ind[$first][$pos];
				$ind[$first][$pos]=$ind[$second][$pos];
				$ind[$second][$pos]=$aux;
			}
		}
	}
#	print "\n";
}

# Mutation
sub mutation {

	my $myrand;

	# For each subject...
	for my $row (0..scalar(@ind)-1){
#		print "row: " . $row . "\n"; 
		# For each column...
		for my $col (0..($indsize-1)){
#			print "col: " . $col . "\n"; 
			$myrand=rand();
#			print "myrand: " . $myrand . "\n"; 
			if ( $myrand <= $mutationprob ) {
				print "Mutation in individual $row at position $col\n";
#				print "before: " . $ind[$row][$col] . "\n";
#				print "cell + 1: " . ($ind[$row][$col] + 1) . "\n";
				$ind[$row][$col]=(($ind[$row][$col] + 1) & 1);
#				print "after: " . $ind[$row][$col] . "\n";
			}
#			print "\n";
		}
	}
	print "\n";
}


# Cross Validation
sub cross {

	my $chosenindnum=0;
	my $percentage=0;
	my $crosspercentage=0;
	my $reporting_crossresultssize=0;
	my $thingy=0;
	my @spotons;
	
	$reporting_crossresultssize = $crossresultssize + 1;

	print "\n>>> Cross validation results <<<\n";

	# The Chosen One
	print "Best individual: ";
	for my $col (0..($indsize-1)){
		print "$ind[scalar(@ind)-1][$col] ";
        }
	print "\n";

	# Calculating the sum of the columns of the best individual, ie., number of technical indicators chosen
	for my $col (0..($indsize-1)){
		$chosenindnum=$chosenindnum+$ind[scalar(@ind)-1][$col];
	}
	print "# of chosen indicators: $chosenindnum\n";

	for my $col (0..($chosenindnum)){
#		print "col: $col\n";
		$spotons[$col] = 0;
#		print "spotons[$col]: $spotons[$col]\n";
	}

	# Counts how many results are reproduced by using the best individual
	for my $row (0..$crossresultssize) {
#		print "row: $row\n";
		my $rights=0;

		for my $col (0..($indsize-1)){
#			print "col: $col\n";
			if ( $ind[scalar(@ind)-1][$col] ){
				if ( $crossresults[$row][$col] == $crossresults[$row][$indsize] ){
					$rights=$rights+1;
				}
			}
		}
#		print "rights: $rights\n";
		for my $col (1..$rights){
			$spotons[$col] = $spotons[$col] + 1;
		}
	}

	for my $numspoton (1..($chosenindnum)){
#		print "numspoton: $numspoton\n";
#		print "crossresults scalar: $crossresultssize\n";
		$crosspercentage = $spotons[$numspoton] / $reporting_crossresultssize;
		$thingy = $crosspercentage;
		$crosspercentage = $crosspercentage * 100;
		$percentage = ( $numspoton / $chosenindnum ) * 100;
		$thingy = $thingy * ( $numspoton / $chosenindnum ) * 100; 
		printf "Times at least $numspoton (%.2f%",$percentage;
		print ") of the chosen indicators predicted correctly $spotons[$numspoton]/$reporting_crossresultssize (";
		printf "%.2f%",$crosspercentage;
		print ") of the crossvalidation cases.";
		if ( $thingyflag ){
			printf " Thingy: %.2f", $thingy;
		}
		print "\n";
	}
	print "\n";
}


################################################
##### MAIN ########## MAIN ########## MAIN #####
################################################

&popresults;
print ">>> Input file loaded\n";
&printarray(@loadedresults);

&holdout;
print ">>> Holdout executed\n\n";
&printgaconditions;
print ">>> Training data\n";
&printarray(@results);
print ">>> Validation data\n";
&printarray(@crossresults);

print ">>> Initializing individuals...\n";
&initialize();

print ">>> First Population\n";
&printarray(@ind);

print ">>> Calc Fitness\n";
$maxval=&calcfitness();

print ">>> Population after fitness\n";
&printarray(@ind);

print ">>> Cross validation with best individual\n";
&cross;

for my $iteration (2..$iterationnumber){

	print "#################################\n";
	print "### Iteration #" . $iteration . "\n";

	print ">>> Selection\n";
	&selection($maxval);

#	print "New population\n";
#	&printarray(@ind);

	print ">>> Crossover\n";
	&crossover;

	print ">>> Mutation\n";
	&mutation;

#	print "New population\n";
#	&printarray(@ind);

	print ">>> Calc Fitness\n";
	$maxval=&calcfitness();

	print ">>> New population after fitness\n";
	&printarray(@ind);

	print ">>> Cross validation with best individual\n";
	&cross;

}

$end = time();

printf("Elapsed time: %.2f\n", $end - $start);

print "Exiting...\n";
print "GA - END\n";

# The answer is 42

exit 0;
