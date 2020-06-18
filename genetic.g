############################################################################
#                            CONSTANTS
############################################################################
YAGS_TEMP.GENETICS:=rec();
YAGS_TEMP.GENETICS.InfoClass := NewInfoClass("YAGSCoilsInfoClass");
##############################




################ CORE TYPES FOR GENETICS #############

NewIndividual := function(obj) #this should be an abstract class
  return rec(obj:=obj);
end;


########## UTILS FOR GENETICS ############

JoinPopulation := function(population1, population2)
   local population,i;

   population := [];

   for i in [1..Size(population1)] do
     Add(population, population1[i]);
   od;

   for i in [1..Size(population2)] do
     Add(population, population2[i]);
   od;

   return population;
end;

SortPopulation := function(population)
  SortBy(population, function(individual)
     return individual.fitness;
  end);
end;

FittestFromPopulation := function(population)
   local fittestIndividual, individual;

   fittestIndividual := population[1];
   for individual in population do
      if individual.fitness>fittestIndividual.fitness then
          fittestIndividual:= individual;
      fi;
   od;
   return fittestIndividual;
end;

IndexOfWeakestFromPopulation := function(population)
  local i,weakestPos, weakestFitness;

  weakestPos := 1;
  weakestFitness := population[weakestPos].fitness;
  for i in [1..Length(population)] do
     if population[i].fitness<weakestFitness then
         weakestFitness:=population[i].fitness;
         weakestPos:=i;
     fi;
  od;
  return weakestPos;
end;


########### CORE GENETICS ##########


ChooseRandomIndividual:=function(useInverseRoulette, population)
   local selectedIndex, fitnessList;

   fitnessList := List(population, ind->ind.fitness);

   if useInverseRoulette=true then
     selectedIndex := ChooseByInverseRoulette(fitnessList);
   else
     selectedIndex := ChooseByRoulette(fitnessList);
   fi;

   return population[selectedIndex];
end;

RandomPopulation := function(populationSize, newRandonIndividualFunc)
  local i, population, individual;

  population := [];
  for i in [1..populationSize] do
    Info(YAGS_TEMP.GENETICS.InfoClass,1, "total of random individuals created:", i, "\r");
    individual:=newRandonIndividualFunc();
    Add(population, individual);
  od;
  Info(YAGS_TEMP.GENETICS.InfoClass,1, "\n");
  return population;
end;

SumFitness := function(oldPopulation)
   local sumFitness, individual;

   sumFitness := 0;

   for individual in oldPopulation do
      sumFitness := sumFitness + individual.fitness;
   od;
   return sumFitness;
end;


CrossCouple := function(crossProb, population, useInverseRoulette, crossFunction)
    local indA, indB, childs;

    indA := ChooseRandomIndividual(useInverseRoulette, population);
    indB := ChooseRandomIndividual(useInverseRoulette, population);

    if FlipCoin(crossProb)=true then
      childs := crossFunction(indA, indB);
    else
      childs := [indA, indB];
    fi;

    return childs;
end;

MutateChilds := function(mutationProb, childs, mutateFunction)
 local child, mutatedChilds;

 mutatedChilds := [];
 for child in childs do
   Add(mutatedChilds, mutateFunction(child));
 od;

 return mutatedChilds;
end;

PredateWeakestFromSubset := function(subsetSize, population,childs)
   local child, individualsIndexes, individualsFitnesses, i, weakestIndividualIndex;

   for child in childs do

	   individualsIndexes := RandomUnrepeatedIntegerSequence(subsetSize,1,Size(population));
	   individualsFitnesses := [];

	   Info(YAGS_TEMP.GENETICS.InfoClass,3, "chosen individuals indexes:", individualsIndexes, "\n");

	   for i in individualsIndexes do
	      Add(individualsFitnesses, population[i].fitness);
	   od;

	   SortParallel(individualsFitnesses, individualsIndexes);

	   weakestIndividualIndex := individualsIndexes[1];
     if population[weakestIndividualIndex].fitness > child.fitness then
	      Info(YAGS_TEMP.GENETICS.InfoClass,3, "will predate individual with index:", weakestIndividualIndex, "\n");
	      population[weakestIndividualIndex] := child;
     fi;
   od;

end;

InitPopulation := function(maxPopulationSize, initialPopulation, newRandonIndividualFunc)
  return JoinPopulation(
                    initialPopulation,
                    RandomPopulation(
                        maxPopulationSize-Size(initialPopulation),
                        newRandonIndividualFunc
                    )
                );
end;

Genetics := function(parameters)
  local params, genNum, childs, population, coupleCount;

  params := ShallowCopy(parameters);
  if not IsBound(params.civilizationId) then
    ## default civilization id
    params.civilizationId:=1;
  fi;

  population := InitPopulation(params.maxPopulationSize, params.initialPopulation, params.newRandonIndividualFunc);

  for genNum in [1 .. params.maxGenerations] do

    for coupleCount in [1 .. params.maxCouples] do

      Info(YAGS_TEMP.GENETICS.InfoClass,3, "crossing couple \n");
	    childs := CrossCouple(params.crossProb, population, params.useInverseRoulette, params.crossFunction);

	    Info(YAGS_TEMP.GENETICS.InfoClass,3, "mutating childs:", childs, " \n");
	    childs := MutateChilds(params.mutationProb, childs, params.mutateFunction);

	    #choose where to put the new childs in the next generation
	    params.predatoryFunction(population, childs);
    od;
    params.onOneGenerationProcessed(genNum, population);
  od;
  params.onAllGenerationsProcessed(population);
  return population;
end;



GeneticsOnCivilizations := function(params)
  local i, augmentsIndex, civilizationNumber, population, augments, weakAugmentIndex, populations, fittest;

  params := ShallowCopy(params);
  augmentsIndex:=0;
  if params.allowEugenics then
    params.numberOfCivilizations := params.numberOfCivilizations + 1;
    augmentsIndex:= params.numberOfCivilizations; #Khaaaaaaan!!! The wrath of Khan is comming!
    Info(YAGS_TEMP.GENETICS.InfoClass,1,"eugenics is enabled, an additional civilization will be added as number:", augmentsIndex, "\n\n");
  fi;

  if not IsBound(params.initialPopulations) then
    params.initialPopulations := List([1..params.numberOfCivilizations],i->[]);
  fi;

  civilizationNumber:=0;
  populations := List(params.initialPopulations, function(initialPopulation)
     local population;

     civilizationNumber:=civilizationNumber+1;
     Info(YAGS_TEMP.GENETICS.InfoClass,1, "initialazing population for civilization:", civilizationNumber, "\n");
     population:=InitPopulation(params.maxPopulationSize, initialPopulation, params.newRandonIndividualFunc);
     Info(YAGS_TEMP.GENETICS.InfoClass,1, "\n");

     return population;
  end);

  while true do

     for civilizationNumber in [1 .. params.numberOfCivilizations] do

        if civilizationNumber=augmentsIndex then
             Info(YAGS_TEMP.GENETICS.InfoClass,1,"\n\n************** civilization:", civilizationNumber, " (augments) ******\n");
        else
             Info(YAGS_TEMP.GENETICS.InfoClass,1,"\n\n************** civilization:", civilizationNumber, " *****************\n");
        fi;
        population:=populations[civilizationNumber];
        params.initialPopulation:=population;
        params.onOneGenerationProcessed:=function(generationNumber, population)
            params.reportAfterOneGenerationProcessedFunction(civilizationNumber, generationNumber, populations);
        end;
        params.onAllGenerationsProcessed:=function(population)
            params.reportAfterAllGenerationsProcessedFunction(civilizationNumber, populations);
        end;
        #evolve population
        population := Genetics(params);

        if params.allowEugenics and civilizationNumber<>augmentsIndex then #will this start "the eugenics war"?? :-)
          #get the fittest from the current population
          fittest := FittestFromPopulation(population);
          augments:=populations[augmentsIndex];
          weakAugmentIndex:=IndexOfWeakestFromPopulation(augments);
          #replace a weak augment with a strong individual for the augments
          if fittest.fitness > augments[weakAugmentIndex].fitness then
            augments[weakAugmentIndex]:=fittest;
          fi;
        fi;

        #@TODO: add some kind of migration if the fitness for some population hasn't changed

     od;

     params.reportAfterAllCivilizationsProcessed(populations);
  od;
end;
