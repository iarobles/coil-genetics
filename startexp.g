LoadPackage("yags");

############################################################################
#                            CONSTANTS
############################################################################

#DeclareGlobalVariable("YAGS_TEMP");
#InstallGlobalVariable(YAGS_TEMP, rec());
YAGS_TEMP := rec();
Read("core.g");
Read("regressions.g");
Read("graph-utils.g");
Read("random-utils.g");
Read("coils-old.g");
Read("genetic.g");
Read("genetic-coils.g");

##########################################################################
########################## GLOBAL VARIABLES ##############################
##########################################################################


  #######################################################################
  ###                 genetic parameters                               ##
  #######################################################################
  #number of simultaneous civilizations being evolved
  NUMBER_OF_CIVILIZATIONS := 3;
  # if true, an extra civilization will be added that will include the
  # most fit individuals of each civilization
  ALLOW_EUGENICS := true;
  MAX_POPULATION_SIZE := 50;
  MAX_GENERATIONS := 20;
  CROSS_PROBABILITY := 1/10;
  MUTATION_PROBABILITY := 3/10;
  MAX_COUPLES := 1; #number of chosen couples that will bear a child into this wild world :p
  #MUTATION PROABILITIES FOR NUMBERS ON INTERVALS
  MUTATE_INTERVAL_ADD_NEW_NUMBER_PROBABILITY := 1/1000;
  MUTATE_INTERVAL_KEEP_NUMBER_PROBABILITY :=  5/10;

  #######################################################################
  ###  parameters for random coil graphs (new individuals)             ##
  #######################################################################
  MAX_RANDOM_N := 30; #should be large enough
  MIN_RANDOM_N := 25;
  MIN_RANDOM_ORBITS := 1;
  MAX_RANDOM_ORBITS := 1;
  #if T^{ab} are the connections from orbit a to orbit b, the following
  #constants ensure the restriction:
  # T^{ab} \subseteq \{MIN_RANDOM_INTERVAL_VAL,...,MAX_RANDOM_INTERVAL_VAL \}
  MIN_RANDOM_INTERVAL_VAL := -8;
  MAX_RANDOM_INTERVAL_VAL := 8;

  #######################################################################
  ###  parameters for fitness                               ##
  #######################################################################
  ## number of iterated clique graphs that will be computed for an individual
  ## in order to compute its fitness
  ## WARNING: greater values produce very slow times
  MAX_ITERATED_CLIQUE_GRAPHS_FOR_FITNESS := 7;
  ## if H is a graph, MAX_CLIQUE_GROW_RATIO ensures the condition:
  ##           |H| < MAX_CLIQUE_GROW_RATIO*|K(H)|.
  ## this restriction is applied for every individual and for each of
  ## computation of its next iterated clique graph, also stops clique computations
  ## before reaching MAX_ITERATED_CLIQUE_GRAPHS_FOR_FITNESS if the condition is
  ## no longer met
  ## WARNING: greater values produce very slow times
  MAX_CLIQUE_GROW_RATIO := 201/100;

  #######################################################################
  ###  name of the report file
  #######################################################################
  FILE_POPULATIONS_REPORT:="reports/populations.rt";
  #if you call Read("reports/populations.rt") it will load the
  #variable "initialPopulationsFromFile" wich contains all the populations
  #that were saved from your previous run of this file

  FILE_POPULATIONS_CONNECTION_MATRIX_REPORT:="reports/populations_connection_matrix.rt";

  FILE_FITTEST_OF_EACH_CIVILIZATION:="reports/fittest_of_each_civilization.rt";

  #######################################################################
  ###                  OPTIONS FOR LOGGING                        #######
  #######################################################################
  SetInfoLevel(YAGS_TEMP.GENETICS.InfoClass, 1);
  #override implementation for Info to avoid print "\n" at the end
  SetInfoHandler(YAGS_TEMP.GENETICS.InfoClass, function(infoclass, level, list)
     local s;
     for s in list do
       Print(s);
     od;
  end);



##########################################################################
#############  FUNCTIONS FOR DEFINING EXPERIMENT #########################
##########################################################################

  #two args are needed due to NumberOfIteratedCliquesUsingRestriction on main.g
  ComputeCliqueGrowthRestriction := function(G,k)
    local growth;

    growth:=Int(Ceil(Float(MAX_CLIQUE_GROW_RATIO)*Order(G)));
    Info(YAGS_TEMP.GENETICS.InfoClass, 5, "for k:",k," computed growth using |G|=",Order(G), " is:", growth,"\n");
    return growth;
  end;


  ############# CORE TYPE for defining a new individual ###################
  NewIndividual := function(coil) #should this function be some kind of abstract class?
     local fitness;
     Info(YAGS_TEMP.GENETICS.InfoClass, 2, "will compute fitness for new coil with connection matrix:", coil.connectionMatrix,"\n");
     fitness := FitnessForIteratedCliquesByLogRegressionUsingRestriction(
                  MAX_ITERATED_CLIQUE_GRAPHS_FOR_FITNESS,
                  coil.G,
                  ComputeCliqueGrowthRestriction
               );
     Info(YAGS_TEMP.GENETICS.InfoClass, 2, "fitness is ", fitness, " for new coil with connection matrix:", coil.connectionMatrix,"\n");
     #param obj and fitness are required to ensure compatibility with functions on genetic.g
     return rec(obj:=coil, fitness:=fitness);
  end;

  SaveExperiment := function(populations)
    local fittestOfEachCivilization, fittest, matrices;

    fittestOfEachCivilization:=[];
    populations:=List(populations,function(population)
       #sort all populations from stronger to weakest before saving
       SortPopulation(population); #@TODO add a flag to reverse the order
       population:=Reversed(population);
       fittest:=population[1];
       Add(fittestOfEachCivilization,rec(connectionMatrix:=fittest.obj.connectionMatrix, fitness:=fittest.fitness));
       return population;
    end);
    matrices:=List(populations,population->List(population,individual->rec(connectionMatrix:= individual.obj.connectionMatrix,fitness:=individual.fitness)));

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "\n\n-------------------------------------------------------\n");

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "saving connection matrix of fittest of each civilization in file:", FILE_FITTEST_OF_EACH_CIVILIZATION,"\n");
    PrintTo(FILE_FITTEST_OF_EACH_CIVILIZATION, "fittestOfEachCivilization:=", fittestOfEachCivilization,";");

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "saving connection matrix of populations in file:", FILE_POPULATIONS_CONNECTION_MATRIX_REPORT,"\n");
    PrintTo(FILE_POPULATIONS_CONNECTION_MATRIX_REPORT,matrices);

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "saving populations in file:", FILE_POPULATIONS_REPORT,"\n");
    PrintTo(FILE_POPULATIONS_REPORT, "initialPopulationsFromFile:=", populations,";");

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "-------------------------------------------------------\n\n");
  end;

##########################################################################
########################## EXECUTION OF EXPERIMENT #######################
##########################################################################

Print("starting experiment using the following parameters:\n\n");
Print("------------------------- core genetic params -------------------------\n");
Print("max number of generations: ",MAX_GENERATIONS,"\n");
Print("population size: ",MAX_POPULATION_SIZE,"\n");
Print("couples: ",MAX_COUPLES,"\n");
Print("cross probability: ",CROSS_PROBABILITY,"\n");
Print("mutation probability: ",MUTATION_PROBABILITY,"\n");

Print("\n----- random coils params (for generating random individuals) -------\n");
Print("number of orbits will be in interval:            [", MIN_RANDOM_ORBITS, ",..., ", MAX_RANDOM_ORBITS, "]\n");
Print("number of vertices per orbit will be in interval:[", MIN_RANDOM_N, ",..., ", MAX_RANDOM_N, "]\n");
Print("values of connection matrix will be in interval: [", MIN_RANDOM_INTERVAL_VAL, ",..., ", MAX_RANDOM_INTERVAL_VAL, "]\n");

Print("\n----------------- fitness params (for each individual) ----------------\n");
Print("maximum number of computed iterated clique graphs:",MAX_ITERATED_CLIQUE_GRAPHS_FOR_FITNESS ,"\n");
Print("maximum clique growth ratio:",MAX_CLIQUE_GROW_RATIO ,"\n\n");

#create params for experiment
geneticsParams := rec(
  numberOfCivilizations:=NUMBER_OF_CIVILIZATIONS,
  allowEugenics:= ALLOW_EUGENICS,
  #initialPopulations := List([1..NUMBER_OF_CIVILIZATIONS],i->[]),
  useInverseRoulette := false,
  maxGenerations := MAX_GENERATIONS,
  maxPopulationSize := MAX_POPULATION_SIZE,
  crossProb := CROSS_PROBABILITY,
  mutationProb := MUTATION_PROBABILITY,
  maxCouples := MAX_COUPLES,
  newRandonIndividualFunc := function()
    local coil;
    coil:=FullRandomCoilGraph(MIN_RANDOM_N,
                        MAX_RANDOM_N,
                        MIN_RANDOM_ORBITS,
                        MAX_RANDOM_ORBITS,
                        MIN_RANDOM_INTERVAL_VAL,
                        MAX_RANDOM_INTERVAL_VAL);
    return NewIndividual(coil);
  end,
  crossFunction := function(indA, indB)
    return [NewIndividual(CrossCoilGraphsReturnOneChild(indA.obj, indB.obj))];
    #return [NewIndividual(CrossCoilGraphs(indA.obj, indB.obj))];
  end,
  mutateFunction := function(individual)
     return NewIndividual(
               MutateCoilGraph(
                  individual.obj,
                  MIN_RANDOM_INTERVAL_VAL,
                  MAX_RANDOM_INTERVAL_VAL,
                  MUTATE_INTERVAL_ADD_NEW_NUMBER_PROBABILITY,
                  MUTATE_INTERVAL_KEEP_NUMBER_PROBABILITY
               )
            );
  end,
  predatoryFunction := function(population, childs)
    #the following function is defined on genetic.g
    PredateWeakestFromSubset(2, population, childs);
  end,
  reportAfterOneGenerationProcessedFunction := function(civilizationNumber, generationNumber, populations)
    #stuff that you may want to report after one generation of civilizationNumber has been processed
    local fittestIndividual;
    #get the fittest individual to report the maximum fitness of the population for this generation
    fittestIndividual := FittestFromPopulation(populations[civilizationNumber]);

    Info(YAGS_TEMP.GENETICS.InfoClass, 1, "generations, total/current: ",MAX_GENERATIONS,"/", generationNumber," max fitness:", Float(fittestIndividual.fitness), "                        \r");
  end,
  reportAfterAllGenerationsProcessedFunction := function(civilizationNumber, populations)
    #stuff that you may want to report after ALL generations for civilizationNumber has been processed
    #SaveExperiment(populations);
  end,
  reportAfterAllCivilizationsProcessed := function(populations)
    #stuff that you may want to report after ALL generations of ALL civilizations has been processed
    SaveExperiment(populations);
  end
);


GeneticsOnCivilizations(geneticsParams);
