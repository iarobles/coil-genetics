
################ COIL GRAPH UTILS FOR NEIGHBORHOOD OF VERTICES  #############

# x should be in [-n/2 .. 0 .. n/2]
# orbnum should be in [1..totorb]
# returns vertex in [1.. totOrb*n];
Cvertex := function(g,orbnum,x)
  local v;

  if x >= 0 then
    v := x + 1;
  else
    v := g!.n + x + 1;
  fi;

  v:= g!.orbits[orbnum][v];

  return v;
end;

# x should be in [1.. totOrb*n]
# returns vertex in [n/2 .. 0.. n/2]
#
Cvertex2 := function(g, orbnum, x)
   local v;

   v := Position(g!.orbits[orbnum],x);

   #change to format [-n/2 .. 0 .. n/2]
   if v>(g!.n)/2 then
     v := v-1-(g!.n);
   else
     v := v-1;
   fi;

   return v;
end;

ToRows:=function(g,q)
  return List(g!.orbits, o->
       Set(List(Intersection(q,o),
           function(c)
              if c-o[1] < g!.n/2 then
                return c-o[1];
              else
                return c-o[1]-g!.n;
              fi;
           end
       ))
   );
end;

#
# x should be in [-n/2 .. 0 .. n/2]
#
CA := function(g, orbnum, x)
  local v, adjs;

  v := Cvertex(g, orbnum, x);
  adjs := Concatenation(Adjacencies(g)[v],[v]);
  adjs := List(g!.orbits, orbit->Intersection(orbit, adjs));

  #PrintListOneItemPerLine(adjs);

  adjs := List([1..Size(g!.orbits)], orbindx->Set(List(adjs[orbindx],v->Cvertex2(g, orbindx, v))));
  #Perform(adjs, Sort);

  return adjs;
end;

CI := function(col)
   local totorbs, colsize, adjs, i;

   totorbs := Size(col[1]);
   colsize := Size(col);

   adjs := [];
   for i in [1..totorbs] do
      Add(adjs, Intersection(List([1..colsize],j-> col[j][i])));
   od;

   return adjs;
end;
