#!/usr/bin/perl

use strict;
use Time::HiRes qw( time );
use POSIX qw/ceil floor/;
use threads;
use threads::shared;
use Thread qw(async);
use List::Util qw[min max sum];
use Math::BigInt;

my $start;
$start = time();

print "GA - START\n";

# Parameters
my $iterationnumber=150;		# number of iterations
my $crossoverprob=0.3;			# crossover probability
my $mutationprob=0.01;			# mutation probability
my $fitnesspow=2;			# pow value for fitness formula
my $indnumber=25;			# number of individuals (solutions)
my $holdoutval=0.1;			# percentage of inputs that will be used for cross validation
my $infile="./input.txt";		# indicators results input file
my $avgenvelope=0.3333;			# prediction_avg envelope
my $initconst=0.15;			# initialize: chance of bit=1

# Variables
my $linecount=-1;
my $maxval=0;
my $end;
my $indsize;			# Size of individuals - number of genes
my $crossresultssize=0;
my $trainingsize=0;
my $loadedresultssize=0;
my $currentiteration=1;		# current iteration
my $onethird=1/3;
my $minusonethird=-$onethird;
my @ind;
my @training;
my @fields;
my @crossresults;
my @crossval;
my @loadedresults;
my $MAX_THREADS = 4;
my $DEBUG_MODE = 0;

my $idxFitness;
my $idxFitnessAdj;
my $idxUsedind;

#system(clear);

# Reading parameters
# Reading input file
if ( $ARGV[0] ){
	$infile = $ARGV[0];
}
# Reading population size. If empty, default value, $indnumber, will be used
if ( $ARGV[1] ){
	$indnumber = $ARGV[1];
}

# Reading number of iterations. If empty, default value, $iterationnumber, will be used
if ( $ARGV[2] ){
	$iterationnumber = $ARGV[2];
}

# Reading avg envelope. If empty, default value, $avgenvelope, will be used
if ( $ARGV[3] ){
	$avgenvelope = $ARGV[3];
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
	print "Number of genes:\t\t$indsize\n";
	print "Input file:\t\t\t$infile\n";
	print "Total data set size:\t\t$loadedresultssize\n";
	print "Training data set size:\t\t$trainingsize\n";
	print "Cross validation data set size:\t$crossresultssize\n";
	print "Avg envelope:\t\t\t$avgenvelope\n";
	print "Init const:\t\t\t$initconst\n";
	print "\n";
}

sub printholdoutconditions {
	print "Training data set size:\t\t$trainingsize\n";
		print "Cross validation data set size:\t$crossresultssize\n";
}

# populates a matrix with the indicator results calculated in the spreadsheet
sub popresults {
	my $line = 0;

	open(INFILE,$infile) or die("Cant open file:$!");
	while (<INFILE>){
		my $line=$_;
		chomp($line);
		@fields = split ',', $line;
		$linecount=$linecount+1;
			for my $colcount (0..scalar(@fields)-1){
			$loadedresults[$linecount][$colcount] = @fields[$colcount];
		}
	}
	close(INFILE);

	# Defining individual size based on input file
	$indsize = scalar(@fields)-1;

	$idxFitness    = $indsize;
	$idxFitnessAdj = $indsize+1;
	$idxUsedind    = $indsize+2;

	# Defining loaded results size
	$loadedresultssize = scalar(@loadedresults);
}

# Holdout method to divide the loadedresults matrix in two parts
sub holdout {
	my $trainingrow=0;
	my $crossvalrow=0;
	@training = ();
	@crossresults = ();
	for my $row (0..scalar(@loadedresults)-1){
		if ( rand() > $holdoutval ){
			$training[$trainingrow]=$loadedresults[$row];
			$trainingrow=$trainingrow+1;
		}
	}
	for my $row (0..scalar(@loadedresults)-1){
		if ( rand() > $holdoutval ){
			$crossresults[$crossvalrow]=$loadedresults[$row];
			$crossvalrow=$crossvalrow+1;
		}
	}
	# Defining training data set size
	$trainingsize = scalar(@training);
	# Defining cross validation data set size
	$crossresultssize = scalar(@crossresults);
}

# Calculates performance of each indicator if used separately
sub indicatorstest {
	my $crossrights=0;
	my $trainingrights=0;
	my $crossratio=0;
	my $trainingratio=0;

	print "Correct predictions per indicator\n";
	print "Indicator\tTraining ($trainingsize)\t\tValidation ($crossresultssize)\n";
	# For each indicator...
	for my $indicator (0..($indsize-1)){
		$crossrights=0;
		$trainingrights=0;
		# For each training validation data set row...
		for my $trainingrow (0..($trainingsize-1)){
			# Compare if it predicted correctly
			if ( $training[$trainingrow][-1] == $training[$trainingrow][$indicator] ){
				$trainingrights++;
			}
		}
		# For each cross validation data set row...
		for my $crossrow (0..($crossresultssize-1)){
			# Compare if it predicted correctly
			if ( $crossresults[$crossrow][-1] == $crossresults[$crossrow][$indicator] ){
				$crossrights++;
			}
		}
		$trainingratio = $trainingrights / $trainingsize;
		$trainingratio = $trainingratio * 100;
		$crossratio = $crossrights / $crossresultssize;
		$crossratio = $crossratio * 100;
		printf "\t$indicator\t$trainingrights\t%.2f%\t\t$crossrights\t%.2f%\n",$trainingratio,$crossratio;
	}
	print "\n";
}

# Print arrays
sub printarray {

	my @array=@_;
	my $sum;

	for my $indn (0..scalar(@array)-1){
		$sum=0.0;
		print "\tLine\t$indn:\t";
#		for my $column (0..scalar(@{$array[0]})-1){
		for my $column (0..$indsize-1){
			printf "%s",$array[$indn][$column];
		}
		printf " | %5s", $array[$indn][$idxFitness];
		printf " | %5s", $array[$indn][$idxFitnessAdj];
		printf " | %2s", $array[$indn][$idxUsedind];
		print "\n";
	}
	print "\n";
}
# Initialize subjects
sub initialize {
	print ">>> Initializing individuals\n";

	my $value;
	my $const = $initconst;

	for my $indn (0..$indnumber-1){
		for my $column (0..$indsize-1){

			my $myrand=rand();
			if ( $myrand gt $const ){
				$value=1;
			}
			else {
				$value=0;
			}
			$ind[$indn][$column]=$value;
		}
		$ind[$indn][$idxFitness]=0.00;
		$ind[$indn][$indsize+1]=0.00;
		$ind[$indn][$indsize+2]=0.00;

		$const = 1 - $const;
	}
}

# Normalize fitness column for next generation selection through roulette wheel method
sub fitnessadj {

	my $max=0.00;
	my $diffsum=0.00;
	my $fitnessval=0.00;
	my $diffabs=0.00;

	# Getting the max of fitness column
	for my $row (0..scalar(@ind)-1){
		$fitnessval=$ind[$row][$idxFitness];
		if ( $max < $fitnessval ) {
			$max=$fitnessval;
		}
	}
	# abs of diff between fitness and its max
	for my $row (0..scalar(@ind)-1){
		$diffabs = abs($ind[$row][$idxFitness]-$max);
		$ind[$row][$idxFitnessAdj] = $diffabs ** $fitnesspow;;
		$ind[$row][$idxUsedind] = usedindnum($row);
		# Getting the sum of (fitness - max)^fitnesspow
		$diffsum += $ind[$row][$idxFitnessAdj];
	}

	return $diffsum;
}

# Returns number of used indicators by individual which row is passed as first parameter
sub usedindnum {
	my $indrow = $_[0];
	return sum(@{$ind[$indrow]}[0..$indsize-1]);
}

# Prediction function for a given indicator and result row. AVG method.
# Returns 2 if number of indicators used equals 0
# Returns 17 in case of errors in the logic
sub prediction_avg {
	my $array = $_[0];
	my $arrayrow = $_[1];
	my $indrow = $_[2];

	my $sum = 0;
	my $avg = 0;
	my $usedind = 0;
	my $pred = 17;

	# Sum of used indicators results
	for my $col (0..($indsize-1)){
		$sum += ( $ind[$indrow][$col] * $array->[$arrayrow][$col] );
	}
	$usedind = usedindnum($indrow);
	if ( $usedind == 0 ) {
		$pred = 2;
	} else {
		$avg = $sum / $usedind;

		# Mutable envelope
		# $avgenvelope = 1 / $usedind;
		$avgenvelope = ceil($usedind / 2.0) / $usedind;

		if (( $avg >= $avgenvelope ) and ( $avg <= 1 )) {
			$pred = 1;
		} else {
			if (( $avg <= ((-1) * $avgenvelope) ) and ( $avg >= -1 )) {
				$pred = -1;
			} else {
				$pred = 0;
			}
		}
	}
	return $pred;
}

sub theFitness_avg {
	my $subject = shift;
	my $sum = 0.00;
	my $myprediction;
	my $error;

	# For each training row...
	for my $trow (0..(scalar(@training)-1)){
		# Calculate absolute error actual result and prediction
		$myprediction = prediction_avg(\@training, $trow, $subject);
		if ( $myprediction == 2 ) {
			# If prediction == 2 (that happens when no indicators are selected)
			# assign worst score possible. "* 2" because diff between 1 and -1 is 2.
			##$sum = $trainingsize * 2;
			$sum = 999999;
			last;
		} else {
			$error = ($training[$trow][-1] != $myprediction);
			#$error = $training[$trow][-1] - $myprediction;
			$sum += abs($error);
			#$sum = $sum * $sum; # Should I use this one?
		}
	}
	return $sum;
}

sub calcfitness_avg1 {
	my $sum;

	# For each subject...
	for my $subject (0..scalar(@ind)-1){
		$sum = theFitness_avg($subject);
		$ind[$subject][$idxFitness] = $sum;
	}

	my $max = &fitnessadj;
	&sortarray();

	return $max;
}

# Fitness function AVG method
sub calcfitness_avg {
	my $error;
	my $max;
	my $myprediction;

	# For each subject...
	##	for my $row (0..scalar(@ind)-1){
	my $limit = scalar(@ind)-1;
	my $chunk = min($limit, 1 + int($limit / $MAX_THREADS));
	my $first = 0;
	my @threadList = ();

	my %theFit :shared;
	%theFit = {};

	while ($first < $limit) {
		my $last = min($first + $chunk, $limit);
		my $thd = async {
			for my $subject ($first..$last) {
				my $sum = theFitness_avg($subject);
				{ lock(%theFit); $theFit{$subject} = $sum; }
			}
		};
		push(@threadList, $thd);
		$first = $last
	}
	map { $_->join(); } @threadList;
	for my $subject (0..$limit){
		unless (exists($theFit{$subject})) {
			print "\t-> INVALID: $subject\n";
			next;
		}
		$ind[$subject][$idxFitness] = $theFit{$subject};
	}
	$max = &fitnessadj;
	&sortarray();
	return $max;
}

#Sort ind array by its last field (adjusted fitness)
sub sortarray {
	my @array = sort cmpfunc @ind;
	@ind=@array;
}
sub cmpfunc {
	# my $size;
	# $size=scalar(@{$ind[0]})-1;
	return (($a->[$idxFitnessAdj] <=> $b->[$idxFitnessAdj]) or
			($b->[$idxUsedind] <=> $a->[$idxUsedind]));
}


# Select new population
sub selection {
	my $max=$_[0];
	my $curfitadj = 0;
	my @indaux = ();

	# loop ($indnumber-2) times to get (same number - 1) subjects for next generation
	for my $count (0..$indnumber-2){
		my $myrand = rand() * $max;
		my $sum = 0;
		# for each subject...
		for my $subject (0..$indnumber-1){
			$curfitadj = $ind[$subject][$idxFitnessAdj];
			if ( ($myrand >= $sum ) and ($myrand <= ($curfitadj + $sum) ) ) {
				@{$indaux[$count]} = @{$ind[$subject]};
				last;  # If found, skip the rest of the loop
			}
			$sum += $curfitadj;
		}
	}
	# Best individual goes automatically to new population
	@{$indaux[$indnumber-1]} = @{$ind[$indnumber-1]};
	@ind = @indaux;
	&sortarray();
}

# Crossover
sub crossover {
	my $pairs=0;
	my $aux;
	my $myrand=0;
	my $first=0;
	my $second=0;
	my $addone=0;
	my $remainder=0;
	$pairs=int($indnumber/2);
	# In case of odd number of individuals, add 50% chances of crosssing over pairs shifted by 1
	# ie., in case of 5 individuals, add 50% chances of crossing over paiasr (2-3) and (4-5)
	# without the code below, only pairs (1-2) and (3-4) would have a chance of crossing over
	$remainder=$indnumber%2;
	if ( $remainder ){
		$myrand=rand();
		if ( $myrand > 0.5 ){
			$addone=1;
		}
	}
	# For possible pairs loop...
	for my $count (1..$pairs) {

		# will crossover happen?
		$myrand=rand();
		if ( $myrand <= $crossoverprob ) {

			$first=(2*$count)-1+$addone;
			$second=(2*$count)-2+$addone;
			my $myrandpos=floor(rand()*($indsize));
			print "\tCrossover between $first and $second at position $myrandpos\n" if $DEBUG_MODE;

			# From position pos till the last data position (scalar-3) swap values
			for my $pos ($myrandpos..($indsize-1)) {
				$aux=$ind[$first][$pos];
				$ind[$first][$pos]=$ind[$second][$pos];
				$ind[$second][$pos]=$aux;
			}
		}
	}
	print "\n" if $DEBUG_MODE;
}

# Mutation
sub mutation {
	for my $subject (0..scalar(@ind)-1){
		for my $col (0..($indsize-1)){
			if ( rand() <= $mutationprob ) {
				print "\tMutation in individual $subject at position $col\n" if $DEBUG_MODE;
				#$ind[$subject][$col]=(($ind[$subject][$col] + 1) & 1);
				$ind[$subject][$col] ^= 1;
			}
		}
	}
	print "\n"  if $DEBUG_MODE;
}

sub indAsBitNum {
	my $str = join('', @_[0..$indsize-1]);
	##print "[$str]\n";
	my $b = Math::BigInt->new('0b' . $str);
	return $b;
}

# Cross Validation
sub crossvalidation {

	my $iteration=$_[0];

	my $chosenindnum=0;
	my $percentage=0;
	my $prediction=0;
	my $crosspercentage=0;
	my $rights=0;
	
	# The Chosen One
	print "I$iteration\tBest sol:\t";
	for my $col (0..($indsize-1)){
		printf "%s", $ind[-1][$col];
	}
	printf " [%s]", &indAsBitNum(@{$ind[-1]});

	# Calculating the sum of the columns of the best individual, ie., number of technical indicators chosen
	$chosenindnum = usedindnum(-1);
	printf "\t# ind: %2d/%2d",$chosenindnum,$indsize;

	# Calculating the sum of the columns of the best individual, ie., number of technical indicators chosen
	my %cardinal = ();
	for my $subject (0..$indnumber-1){
		my @arSubject = @{$ind[$subject]};
		my $indAsNumber = &indAsBitNum(@arSubject);
		$cardinal{$indAsNumber}++;
	}
	printf "\t# dist subj: %2d/%2d",scalar keys %cardinal, $indnumber;

	# Counts how many results are reproduced by the best individual
	# For each cross results record

	for my $crow (0..($crossresultssize-1)) {
		$prediction = prediction_avg(\@crossresults,$crow,-1);
		if ( $prediction == $crossresults[$crow][$idxFitness] ){
			$rights++;
		}
	}

	# for my $crow (0..($trainingsize-1)) {
	# 	$prediction = prediction_avg(\@training, $crow, -1);
	# 	if ( $prediction == $training[$crow][$idxFitness] ){
	# 		$rights++;
	# 	}
	# }

	# for my $crow (0..($loadedresultssize-1)) {
	# 	$prediction = prediction_avg(\@loadedresults, $crow, -1);
	# 	if ( $prediction == $loadedresults[$crow][$idxFitness] ){
	# 		$rights++;
	# 	}
	# }

	$crosspercentage = $rights / $crossresultssize;
	# $crosspercentage = $rights / $loadedresultssize;
	$crosspercentage = $crosspercentage * 100;
	printf "\tCorrect $rights/$crossresultssize (or %.2f%)\n",$crosspercentage;
	# printf "\tCorrect $rights/$loadedresultssize (or %.2f%)\n",$crosspercentage;
	print "\n\n";
}


################################################
##### MAIN ########## MAIN ########## MAIN #####
################################################

&popresults;
print ">>> Input file loaded\n";
#&printarray(@loadedresults);

&holdout;
print ">>> Holdout executed\n\n";
&printgaconditions;

#print ">>> Training data\n";
#&printarray(@training);
#print ">>> Validation data\n";
#&printarray(@crossresults);

print ">>> Performance of each indicator if used separately\n";
&indicatorstest;

&initialize();

print ">>> First Population\n";
&printarray(@ind);

print "##########################################\n";
print "### Iteration #1\n\n";

print "I$currentiteration\tCalc Fitness\n";
$maxval=&calcfitness_avg();

print "I1\tPopulation after fitness\n";
&printarray(@ind);

print "I$currentiteration\tCross validation results\n";
&crossvalidation($currentiteration);

for my $currentiteration (2..$iterationnumber){

	print "##########################################\n";
	print "### Iteration #" . $currentiteration . "\n\n";

	print "I$currentiteration\tSelection\n";
	&selection($maxval);

#	print "New population\n";
#	&printarray(@ind);

	print "I$currentiteration\tCrossover\n";
	&crossover;

	print "I$currentiteration\tMutation\n";
	&mutation;

#	print "New population\n";
#	&printarray(@ind);

	&holdout;
	print ">>> Holdout executed\n\n";
	&printholdoutconditions;

	print "I$currentiteration\tCalc Fitness\n";
	$maxval=&calcfitness_avg();

	print "I$currentiteration\tNew population after fitness\n";
	&printarray(@ind);

	print "I$currentiteration\tCross validation with best individual - results\n";
	&crossvalidation($currentiteration);

}

$end = time();

print "##########################################\n";
printf("Elapsed time: %.2f\n", $end - $start);

print "Exiting...\n";
print "GA - END\n";

warn sprintf "%s\n", &indAsBitNum(@{$ind[-1]});

# The answer is 42

exit 0;
