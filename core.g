
############################################################################
#                            CONSTANTS
############################################################################

#DeclareGlobalVariable("YAGS_TEMP");
#InstallGlobalVariable(YAGS_TEMP, rec());
YAGS_TEMP := rec();

YAGS_TEMP.RANDOM_RESOLUTION := 10000;  #for random utilites


##################  DIRTY PATCHES FOR GAP 4.11 ##############################
if not IsBoundGlobal("PrintListOneItemPerLine") then BindGlobal("PrintListOneItemPerLine",PrintArray); fi;

################## UTILS FOR INTEGERS ##########################
SumOfIntegers:=function(n)
  return n*(n+1)/2;
end;
