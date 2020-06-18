

############################################################################
#                            GRAPH UTILITIES
############################################################################


############################################################################
K:=CliqueGraph;
############################################################################

CliqueGraphK := function(k, G)
  local KG, i;

  KG := G;
  for i in [1..k] do
     KG := CliqueGraph(KG);
  od;

  return KG;
end;

#returns [|G|, |K(G)|, |K^2(G)|, ..., |K^k(G)| ]
NumberOfIteratedCliques := function(k, G)
  local listKG, list, i;

  listKG := [G];
  for i in [1..k-1] do
     Add(listKG, CliqueGraph(listKG[Size(listKG)]));
  od;

  list := List(listKG,Order);

  if k>=1 then
    Add(list,NumberOfCliques(listKG[Size(listKG)]));
  fi;

  return list{[2..Size(list)]};
end;


NumberOfIteratedCliquesUsingRestriction := function(k, G, restrFn)
  local listKG, listTotalCli, forcedToStop, i, H, maxCli, tCliques;

  listTotalCli := [];
  listKG := [G];
  forcedToStop := false;

  for i in [1..k] do
     H:=listKG[Size(listKG)];
     maxCli := restrFn(H,i);
     tCliques := NumberOfCliques(H,maxCli);

     Info(YAGSInfo.InfoClass,9, "tCliques:", tCliques,  ", maxCli:", maxCli, "\n");

     if tCliques = maxCli then #if the computation of cliques was forced to stop
       forcedToStop := true;
     else
       if i<k then
          Add(listKG, CliqueGraph(H));
       fi;
       Add(listTotalCli, tCliques);
     fi;
  od;

  return listTotalCli;
end;


############################################################################
##
#F  ClosedAdjacency ( <G>, <vertex> )
##
##  <Description>
##
##     Returns the closed neighborhood of the vertex <vertex> in the graph
##     <G>, that is:
##        $$ N_{<A>G</A>}[x] = N_{<A>G</A>}(x) \cup \{x\}$$
##
##  </Description>
##
ClosedAdjacency := function(G,x)
    return Union(Adjacency(G,x), [x]);
end;

############################################################################
##
#F  ClosedAdjacencyOfSet ( <G>, <vertexSet> )
##
##  <Description>
##
##     Returns the intersection of the closed neighborhood of each vertex of
##     <vertexSet> in the graph <G>, that is:
##        $$ N_{<A>G</A>}[X] = \bigcap_{x \in X} N_{<A>G</A>}[x]
##     <G>.
##
##  </Description>
##
ClosedAdjacencyOfSet := function(G, vertexSet)
    local NX, x;

    NX := Vertices(G);
    for x in vertexSet do
        NX := Intersection(NX, ClosedAdjacency(G,x));
    od;

    return NX;
end;


############################################################################
##
#F  IsUniversal ( <G> )
##
##  <Description>
##
##     Returns true if the graph <A>G</A> has an universal vertex (i.e
##     a vertex x such that N[x]=G), otherwise returns false.
##
##  </Description>
##
IsUniversal := function(G)
   local top, x;

   top := Order(G) - 1;
   for x in Vertices(G) do
      if Length(Adjacency(G,x)) = top then
         return true;
      fi;
   od;

   return false;
end;

############################################################################
##
#F  IsInducedSubgraph ( <H>, <G> )
##
##  <Description>
##
##     Returns "true" if the graph <H> is an induced subgraph of the graph
##     <G> (i.e H=G[V(H)]).
##
##  </Description>
##
IsInducedSubgraph:=function(h,g)
  return FullMonoMorphism(h,g)<>fail;
end;

############################################################################
##
#F  CircleProduct ( <G>, <H> )
##
##  <Description>
##
##     Returns the circle product of the graphs <G> and <H>.
##
##  </Description>
##
CircleProduct:=function(G,H)
  return ComplementGraph(TimesProduct(ComplementGraph(G), ComplementGraph(H)));
end;


############################################################################
##
#F  FishProduct ( <G>, <H> )
##
##  <Description>
##
##     Returns the fish product of the graphs <G> and <H>.
##
##  </Description>
##
FishProduct := function(G,H)
  local vertices;

  vertices := Cartesian(Vertices(G),Vertices(H));

  return GraphByRelation(vertices, function(x,y)

     if x[1]=y[1] or ( x[1] in Adjacencies(G)[y[1]] and not x[2] in Adjacencies(H)[y[2]] ) then
        return true;
     else
        return false;
     fi;
  end);
end;

############################################################################
##
#F  IsConnectedGraph ( <G>)
##
##  <Description>
##
##      Returns true if the graph <G> is connected and false otherwise.
##
##  </Description>
##
IsConnectedGraph := function(G)
  if NumberOfConnectedComponents(G) = 1 then
    return true;
  else
    return false;
  fi;
end;
