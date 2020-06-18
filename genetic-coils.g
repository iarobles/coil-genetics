

########### UTILS FOR SORTING ##########

SortCoilGraphByOrbitSize := function(coils)
  SortBy(coils, function(coil)
     return Size(coil.connectionMatrix);
  end);
end;


########### UTILS FOR COILS GRAPHS ################

SetConnectionMatrixItem := function(m, row, col, item)
    m[row][col] := item;
    if row<>col then
      m[col][row] := -1*m[row][col];
      Sort(m[col][row]);
    fi;
end;




########### UTILS FOR RANDOM COIL GRAPHS ################

RandomConnectionMatrix := function(orbits, minVal, maxVal)
   local row, col, m, maxRandomNumbers;

   maxRandomNumbers := maxVal - minVal + 1;
   m := NullMat(orbits, orbits);
   for row in [1..orbits] do
     for col in [row..orbits] do
        SetConnectionMatrixItem(m,row,col, RandomUnrepeatedIntegerSequence(maxRandomNumbers, minVal, maxVal));
     od;
   od;
   return m;
end;


RandomCoilGraph := function(verticesPerOrbit, orbits, minVal, maxVal)
  local G, connectionMatrix, coil;
  connectionMatrix := RandomConnectionMatrix(orbits,minVal, maxVal);
  coil:=CoilGraphByConnectionMatrix(verticesPerOrbit, connectionMatrix);
  return coil;
end;

FullRandomCoilGraph :=function(minN, maxN, minOrbits, maxOrbits, minVal, maxVal)
  return RandomCoilGraph(Random(minN,maxN), Random(minOrbits,maxOrbits), minVal, maxVal);
end;


RandomConnectedCoilGraph := function(verticesPerOrbit, orbits, minVal, maxVal)
  local co, connectionMatrix;
  repeat
    connectionMatrix := RandomConnectionMatrix(orbits,minVal, maxVal);
    co := CoilGraphByConnectionMatrix(verticesPerOrbit, connectionMatrix);
  until IsConnectedGraph(co.G)=true;

  return co;
end;

FullRandomConnectedCoilGraph :=function(minN, maxN, minOrbits, maxOrbits, minVal, maxVal)
  return RandomConnectedCoilGraph(Random(minN,maxN), Random(minOrbits,maxOrbits), minVal, maxVal);
end;


############### UTILS FOR CROSSING COIL GRAPHS #############

CrossDiagonalOrbits := function(coilA, coilB)
  local diagonalsToCross, coilATotOrbits, coilBTotOrbits, i, childAm, childBm, diagA, diagB;

  coilATotOrbits := Size(coilA.connectionMatrix);
  coilBTotOrbits := Size(coilB.connectionMatrix);

  childAm := StructuralCopy(coilA.connectionMatrix);
  childBm := StructuralCopy(coilB.connectionMatrix);

  diagonalsToCross := Random(1, coilATotOrbits);
  Info(YAGSInfo.InfoClass,5, "diagonals to cross:", diagonalsToCross, "\n");

  for i in [1..diagonalsToCross] do
    diagA := Random(1,coilATotOrbits);
    diagB := Random(1,coilBTotOrbits);

    Info(YAGSInfo.InfoClass,5, "interchanging (both ways) mA[", diagA, "][", diagA, "]:", coilA.connectionMatrix[diagA][diagA]," with mB[", diagB, "][", diagB, "]", coilB.connectionMatrix[diagB][diagB],"\n");

    SetConnectionMatrixItem(childAm, diagA, diagA, coilB.connectionMatrix[diagB][diagB]);
    SetConnectionMatrixItem(childBm, diagB, diagB, coilA.connectionMatrix[diagA][diagA]);

  od;

  return [CoilGraphByConnectionMatrix(coilA.verticesPerOrbit, childAm), CoilGraphByConnectionMatrix(coilB.verticesPerOrbit, childBm)];
end;

CrossNonDiagonalOrbits := function(coilA, coilB)
  local nonDiagonalsToCross, coilATotOrbits, coilBTotOrbits, i, childAm, childBm, rowA, colA, rowB, colB, childs;

  coilATotOrbits := Size(coilA.connectionMatrix);
  coilBTotOrbits := Size(coilB.connectionMatrix);

  if coilATotOrbits>1 and coilBTotOrbits>1 then

	  childAm := StructuralCopy(coilA.connectionMatrix);
	  childBm := StructuralCopy(coilB.connectionMatrix);

	  nonDiagonalsToCross := Random(1, SumOfIntegers(coilATotOrbits)-coilATotOrbits);
	  Info(YAGSInfo.InfoClass,5, "non diagonals to cross:", nonDiagonalsToCross, "\n");

	  for i in [1..nonDiagonalsToCross] do

	    rowA := Random(1,coilATotOrbits-1);
	    colA := Random(rowA+1,coilATotOrbits);
	    rowB := Random(1,coilBTotOrbits-1);
	    colB := Random(rowB+1,coilBTotOrbits);

	    Info(YAGSInfo.InfoClass,5, "Interchanging (both ways) mA[", rowA, "][", colA, "]:", coilA.connectionMatrix[rowA][colA], " with mB[", rowB, "][", colB, "]:", coilB.connectionMatrix[rowB][colB],"\n");

	    SetConnectionMatrixItem(childAm, rowA, colA, coilB.connectionMatrix[rowB][colB]);
	    SetConnectionMatrixItem(childBm, rowB, colB, coilA.connectionMatrix[rowA][colA]);

	  od;

	  childs := [CoilGraphByConnectionMatrix(coilA.verticesPerOrbit, childAm), CoilGraphByConnectionMatrix(coilB.verticesPerOrbit, childBm)];
  else
      childs := [coilA, coilB];
  fi;

  return childs;
end;


CrossCoilGraphs := function(coilA, coilB)
  local coils, childs;

  coils := [coilA, coilB];
  SortCoilGraphByOrbitSize(coils);
  childs := CrossDiagonalOrbits(coils[1],coils[2]);
  childs := CrossNonDiagonalOrbits(childs[1],childs[2]);

  return childs;
end;

CrossCoilGraphsReturnOneChild := function(coilA, coilB)
  return CrossCoilGraphs(coilA, coilB)[1];
end;


##########################################################################
########################## MUTATE COIL GRAPHS ############################
##########################################################################

MutateIntervalSwapValues := function(interval, minVal, maxVal,addNewNumberProbability, keepNumberProbability)
  local mutatedInterval, number;

  mutatedInterval := [];
  for number in [minVal .. maxVal] do
      if number in interval then
        if FlipCoin(keepNumberProbability)=true  then
           Add(mutatedInterval, number);
        fi;
      else
        if FlipCoin(addNewNumberProbability)=true then
           Add(mutatedInterval,number);
        fi;
      fi;
  od;

  return mutatedInterval;
end;

MutateInterval := function(interval, minVal, maxVal, addNewNumberProbability, keepNumberProbability)
  local mutatedInterval;

  mutatedInterval := MutateIntervalSwapValues(interval, minVal, maxVal, addNewNumberProbability, keepNumberProbability);

  return mutatedInterval;
end;

MutateDiagonalOrbits := function(coil, minVal, maxVal,addNewNumberProbability, keepNumberProbability)
  local m, diagonalsToMutate, orbits, i, diag;

  m := StructuralCopy(coil.connectionMatrix);
  orbits := Size(coil.connectionMatrix);

  diagonalsToMutate := Random(1, orbits);
  Info(YAGSInfo.InfoClass,5, "diagonals to mutate:", diagonalsToMutate, "\n");
  for i in [1..diagonalsToMutate] do
     diag := Random(1,orbits);
     Info(YAGSInfo.InfoClass,5, " mutating m[", diag, "][", diag, "]:=", m[diag][diag], "\n");
     SetConnectionMatrixItem(m,
                             diag,
                             diag,
                             MutateInterval(
                                  m[diag][diag],
                                  minVal,
                                  maxVal,
                                  addNewNumberProbability,
                                  keepNumberProbability
                             )
     );
  od;

  return CoilGraphByConnectionMatrix(coil.verticesPerOrbit, m);
end;

MutateNonDiagonalOrbits := function(coil, minVal, maxVal, addNewNumberProbability, keepNumberProbability)
  local m, nonDiagonalsToMutate, orbits, i, diag, row, col, mutatedCoil, interval;

  orbits := Size(coil.connectionMatrix);

  if orbits > 1 then
      m := StructuralCopy(coil.connectionMatrix);
	  nonDiagonalsToMutate := Random(1, orbits);
	  Info(YAGSInfo.InfoClass,5, "non diagonals to mutate:", nonDiagonalsToMutate, "\n");

	  for i in [1..nonDiagonalsToMutate] do
	     row := Random(1,orbits-1);
	     col := Random(row+1, orbits);

	     interval := MutateInterval(m[row][col], minVal, maxVal, addNewNumberProbability, keepNumberProbability);
	     SetConnectionMatrixItem(m, row, col, interval);

	  od;

	  mutatedCoil := CoilGraphByConnectionMatrix(coil.verticesPerOrbit, m);
  else
      mutatedCoil := coil;
  fi;

  return mutatedCoil;
end;

MutateCoilGraph := function(coil, minVal, maxVal,addNewNumberProbability, keepNumberProbability)
   local mutatedCoil;

   if FlipCoin(1/2) then
     mutatedCoil := MutateDiagonalOrbits(coil, minVal, maxVal, addNewNumberProbability, keepNumberProbability);
     mutatedCoil := MutateNonDiagonalOrbits(mutatedCoil, minVal, maxVal, addNewNumberProbability, keepNumberProbability);
   else
     #Print("\n\nCLIQUE!, n:", coil.verticesPerOrbit, ", r:", coil.totalOrbits, ",m:", coil.connectionMatrix, "\n\n");
     mutatedCoil := CoilCliqueGraph(coil);
   fi;

   return mutatedCoil;
end;


######################## FITNESS FUNCTIONS ###################

FitnessForCliqueGrowRatio := function(G, cliqueGrowRatio, leftOffset, rigthOffset, punishValue)
   local minCliqueGrowRatio, maxCliqueGrowRatio, minCli, maxCli, cliquesCounted, fitness, forcedToStop;

    minCliqueGrowRatio := cliqueGrowRatio - leftOffset;
    maxCliqueGrowRatio := cliqueGrowRatio + rigthOffset;

    minCli := Int(Floor(Float(minCliqueGrowRatio*Order(G))));
    maxCli := Int( Ceil(Float(maxCliqueGrowRatio*Order(G))));

    cliquesCounted := NumberOfCliques(G, maxCli);

    forcedToStop := false;
    if cliquesCounted=maxCli then
      forcedToStop := true;
    fi;

    #option 1 for fitness:
    #fitness := cliqueGrowRatio - AbsoluteValue(cliqueGrowRatio - (cliquesCounted/Order(G)));

    #option 2 for fitness:
    #if  minCli<=cliquesCounted and cliquesCounted < maxCli then
    fitness := cliqueGrowRatio - AbsoluteValue(cliqueGrowRatio - (cliquesCounted/Order(G)));
    if cliquesCounted = maxCli then
       fitness := fitness*punishValue; #punish values that exceed the expected value
    fi;

    Info(YAGSInfo.InfoClass,9,rec(
      minCliqueGrowRatio:=Float(minCliqueGrowRatio),
      maxCliqueGrowRatio:=Float(maxCliqueGrowRatio),
      minCli:=minCli,
      maxCli:=maxCli,
      cliquesCounted:=cliquesCounted,
      fitness:=Float(fitness),
      forcedToStop:=forcedToStop
      ),
      " \n");

    return rec(fitness:=fitness, forcedToStop:=forcedToStop);
end;

FitnessForIteratedCliqueGrowRatio := function(G, maxIterations, cliqueGrowRatio, leftOffset, rigthOffset, punishValue)
   local currentGraph, fitness, fitnessInfo, i;

   currentGraph := G;

   #option 1 for computing fitness (grows exponentially, but there is a problem with forcetToStop)
   fitness := 1;

   #option 2 (grows linearly)
   #fitness := 0;

   for i in [1..maxIterations] do

      Info(YAGSInfo.InfoClass,9,"computing pair , K",i-1,"G, K",i,"G \n");
      fitnessInfo := FitnessForCliqueGrowRatio(currentGraph, cliqueGrowRatio, leftOffset, rigthOffset,punishValue);

      #option 1 for computing fitness (grows exponentially, but there is a problem with forcetToStop)
      fitness := fitness*fitnessInfo.fitness;

      #option 2 (grows linearly)
      #fitness := fitness + fitnessInfo.fitness;

      if fitnessInfo.forcedToStop = false then
         currentGraph := CliqueGraph(currentGraph);
      else
         break; #stop iterations
      fi;

   od;

    Info(YAGSInfo.InfoClass,9, " computed fitness:=", fitness, " \n");

   return fitness;
end;

CountRepeatedAtListEnd := function(list)
   local repeated, i, listSize;

   listSize := Size(list);
   repeated := 0;

   if listSize > 1 then
     for i in [1 .. (listSize-1)] do
       if list[(listSize-i+1)]=list[listSize-i] then
         repeated := repeated + 1;
       else
         break;
       fi;
     od;
   fi;

   if repeated >= 1 then
     repeated := repeated + 1;
   fi;

   return repeated;
end;


FitnessForIteratedCliquesByLogRegressionUsingRestriction := function(maxIterations,G, restrFn)
   local yList, yListSize, xList,r, fitness, regrFitness,nonConstfitness, cliqueIterationFitness;

   yList:= Concatenation([Order(G)], NumberOfIteratedCliquesUsingRestriction(maxIterations,G, restrFn));
   yListSize := Size(yList);

   if yListSize >= 1 then
       xList:= [0 .. (yListSize-1)];

       # 0 = bad aproximation,
       # 1 = perfect aproximation
       if Float(_Variance(xList))=Float(0) or Float(_Variance(yList))=Float(0) then
         regrFitness := Float(1);  # avoids inf problem when computing rSquared
                                   # for xList=[oneValue] (small, only one item) or
                                   # yList=[constant ... constant]
                                   #

       else
         r := LogarithmicRegression(xList,yList);
         #Print("r:",r,", x:", xList, "\n", "y:", yList, "\n");

         if r.b <=Float(1) then  #punish graphs with negative exponential growth
            regrFitness := Float(0);
         else
            regrFitness := r.rSquared;
         fi;
       fi;

       # 0 = yList is constant,
       # 1 = non constant at all
       nonConstfitness := Float(1-(CountRepeatedAtListEnd(yList)/yListSize) );

       # 0 = no iterations were performed,
       # 1 = all iterations were performed
       cliqueIterationFitness := Float((yListSize-1)/maxIterations);
       fitness := regrFitness * cliqueIterationFitness * nonConstfitness;
       Info(YAGS_TEMP.GENETICS.InfoClass,
            2,
            "computed fitness (detail):\n",
            rec(
               yList:=yList,
               fitness:=fitness,
               regrFitness:=regrFitness,
               nonConstfitness:=nonConstfitness,
               cliqueIterationFitness:=cliqueIterationFitness
             ),
            "\n"
       );
   else
       fitness := Float(0);
   fi;

   return fitness;
end;

ConnectedGraphFitness := function(G)
  if IsConnectedGraph(G) = true then
    return 1;
  else
    return 0;
  fi;
end;
