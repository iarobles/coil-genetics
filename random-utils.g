


########### UTILS FOR RANDOM OPERATIONS #####################

RandomUnrepeatedIntegerSequence := function(maxRandomNumbers, minVal, maxVal)
  local sequence;

  sequence := List([1..maxRandomNumbers], i-> Random(minVal, maxVal));
  #Intersection also sorts the sequence
  sequence := Intersection(sequence, sequence);

  return sequence;
end;

RandomRational := function(maxRational)
  return (Random(1, YAGS_TEMP.RANDOM_RESOLUTION)/YAGS_TEMP.RANDOM_RESOLUTION)*maxRational;
end;


FlipCoin := function(probability)

   if RandomRational(1) <= probability then
     return true;
   else;
     return false;
   fi;
end;



############################################################################
##
#F  ChooseByRoulette( <numList> )
##
##  <Description>
##
##     Returns a index from the list <A>numList</A> of a positive
##     real number choosen by using a random roulette (bigger
##     numbers have less probability of been selected).
##
##  </Description>
##
ChooseByRoulette := function(numList)
  local x, sum, rand, partSum, i;

  for x in numList do
    if Float(x)<Float(0) then
        Error("numList is not a list of positive real numbers");
    fi;
  od;

  sum := Sum(numList);
  rand := RandomRational(sum);
  partSum := 0;

  i:=0;
  repeat
    i:=i+1;
    partSum := partSum + numList[i];
  until partSum>=rand or i=Size(numList);

  Info(YAGSInfo.InfoClass,9, "random roulette chose index:", i,  "\n");

   return i;
end;

############################################################################
##
#F  ChooseByInverseRoulette( <numList> )
##
##  <Description>
##
##     Returns a index from the list <A>numList</A> of a positive
##     real number choosen by using a inverse random roulette
##     (bigger numbers have less probability of been selected).
##
##  </Description>
##
ChooseByInverseRoulette := function(numList)

  Info(YAGSInfo.InfoClass,9, "inverse ");
  return ChooseByRoulette(1 - (numList/Sum(numList)) );
end;
