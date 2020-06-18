


############################################################################
#                               COIL GRAPHS
############################################################################



Connections2 := function(C,i,j)
  local conn;

  conn := List(C.connectionMatrix[i][j], function(k)
      local half;

      half := C.verticesPerOrbit/2;
      if k > half then
         return k-C.verticesPerOrbit;
      else
         return k;
      fi;
   end);
   Sort(conn);

   return conn;
end;



CoilGraph := function(record) #@TODO: define a GAP type for this entity

  #@TODO: add some verifications for the following properties in the record:
  #  G:=G,
  #  verticesPerOrbit:=n,
  #  gamma := gamma,
  #  totalOrbits := totalOrbits,
  #  orbits := Orbits(Group(gamma), Vertices(G)),
  #  connectionMatrix := connectionMatrix (integers)

  record.connectionMatrix := List([1..record.totalOrbits], i-> List([1..record.totalOrbits], j-> Connections2(record,i,j)));

  return record;
end;


EmptyConnectionMatrix := function(totalOrbits)
  local m,row,col, vector;

  m:=[];
  for row in [1..totalOrbits] do
    vector := [];
    for col in [1 .. totalOrbits] do
      Add(vector, []);
    od;
    Add(m,vector);
  od;

  return m;
end;

IsSquareMatrix := function(m)
  local row, totalRows, totalColumns;

  totalRows := Size(m);
  for row in [1..totalRows] do
    if Size(m[row])<>totalRows then
       return false;
    fi;
  od;
  return true;
end;

IntegerMatrixToZnMatrix := function(Zn, m)
  local row, col, m2, vector, interval;
  m2 := [];

  for row in [1.. Size(m)] do
    vector := [];
    for col in [1.. Size(m)] do
      #Print("row:", row, ", col:", col, "\n");
      interval := m[row][col]*One(Zn);
      Sort(interval);
      Add(vector, interval);
    od;
    Add(m2,vector);
  od;

  return m2;
end;


ZnMatrixToIntegerMatrix := function(Zn, m)
  local row, col, m2, vector, interval;
  m2 := [];

  for row in [1.. Size(m)] do
    vector := [];
    for col in [1.. Size(m)] do
      interval := List(m[row][col], Int);
      Sort(interval);
      Add(vector, interval);
    od;
    Add(m2,vector);
  od;

  return m2;
end;


############################################################################
##
#F  FormatConnectionMatrix( <m> )
##
##  <Description>
##
##     Returns a new matrix m2 of the same size that m, such that
##          m2[i][j]=m[i][j] for i<>j and
##          m2[i][i]=Union(m[i][j],-m[i][j])
##     This operation ensures that the definition of adjacencies for coil graphs
##     is consistent (i.e. x~y if and only if y~x).
##  </Description>
##
FormatConnectionMatrix := function(m)
   local m2,i,j;

   m2 := [];
   for i in [1..Size(m)] do
     m2[i] := [];
     for j in [1..Size(m)] do
        if i=j then
           m2[i][j]:=Union(m[i][j],-m[i][j]);
        else
           m2[i][j]:=Union([],m[i][j]);
        fi;
     od;
   od;

   return m2;
end;

VerifyCoilConnectionMatrix := function(connectionMatrix)
  local row, column, Tab, minusTba;

  if IsSquareMatrix(connectionMatrix)=false then
    Error("the given connection Matrix is not a square matrix");
  fi;

  for row in [1.. Size(connectionMatrix)] do
     for column in [1.. Size(connectionMatrix)] do

       if row<>column then

           Tab := connectionMatrix[row][column];
           minusTba := -1*Tab;
           Sort(minusTba);

           if connectionMatrix[column][row] <> minusTba then
             Error("\n Matrix is not consistent in row ", row, " and column ", column, ": \n",
             "Matrix[", row, "][", column, "]: ", connectionMatrix[row][column], "\n",
             "Sort of -1*Matrix[", row, "][", column, "]: ", minusTba, "\n",
             "Matrix[", column, "][", row, "]:", connectionMatrix[column][row], "\n");
           fi;
       fi;
     od;
  od;
end;


CoilGraphByConnectionMatrix := function(n, m)
   local vertices, m2, Zn, totalOrbits, interval,gamma,i,G, coil;
   Zn := ZmodnZ(n);
   totalOrbits := Size(m);
   m2:=IntegerMatrixToZnMatrix(Zn, m);
   m2:=FormatConnectionMatrix(m2);
   VerifyCoilConnectionMatrix(m2);

   vertices := Cartesian([0.. (Size(m2)-1)],[0..(n-1)]);
   G:=GraphByRelation(vertices, function(x,y)
      local xOrbit, yOrbit, xIndex, yIndex, distance, adjacent;

      xOrbit := x[1];
      xIndex := x[2];
      yOrbit := y[1];
      yIndex := y[2];
      distance := yIndex*One(Zn) - xIndex*One(Zn);

      if distance in m2[xOrbit+1][yOrbit+1] then
        adjacent := true;
      else
        adjacent := false;
      fi;

      return adjacent;
   end);

  gamma := CycleFromList([1..n]);
  for i in [2..totalOrbits] do
    interval := [((i-1)*n+1)..i*n];
    gamma := gamma*CycleFromList(interval);
  od;

  coil:= CoilGraph(rec(
    G:=G,
    verticesPerOrbit:=n,
    totalOrbits := totalOrbits,
    orbits := Orbits(Group(gamma), Vertices(G)),
    gamma := gamma,
    connectionMatrix := ZnMatrixToIntegerMatrix(Zn, m2)
  ));
  return coil;
end;




############################################################################
##
#F  CirculantCliqueOrbitAdjacencies( <qa>, <qb> )
##
##  <Description>
##
##     Returns the adjacencies (a list) of the clique orbit <A>qa</A> with
##     the clique orbit <A>qb</A> for a circulant graph (i.e. T^{ab}, where
##     T is the connection matrix of the clique graph of a circulant graph).
##
##  </Description>
##
CirculantCliqueOrbitAdjacencies := function(qa,qb)
    local Tab, x, y;

    Tab := [];
    for x in qa do
       for y in qb do
          Add(Tab, x-y);
       od;
    od;

    return Union([],Tab); #cheap trick to order and filter repeated values
end;


############################################################################
##
#F  CoilCliqueOrbitAdjacencies( <qa>, <qb> )
##
##  <Description>
##
##     Returns the adjacencies (a matrix) of the clique orbit <A>qa</A> with
##     the clique orbit <A>qb</A> for a coil graph (i.e. T^{ab}, where
##     T is the connection matrix of the clique graph of a coil graph).
##
##     The clique orbits <A>qa</A> and <A>qb</A> must be matrices (with the
##     same number of rows) such that the $i$-th row is a list of the
##     vertices of the $i$-th orbit of the coil graph that belongs to clique
##     orbit (<A>qa</A> or <A>qb</A>).
##
##  </Description>
##
CoilCliqueOrbitAdjacencies := function(qa,qb)
   local Tab, row, circulantClique;

   Tab := [];
   if Size(qa)<>Size(qb) then
      Error(qa, " and ", qb, " aren't the same size");
   fi;

   for row in [1..Size(qa)] do
      Tab := Union(Tab, CirculantCliqueOrbitAdjacencies(qa[row],qb[row]));
   od;

   return Tab;
end;




######################################################################
##
#F  ConnectionMatrixOfCoilGraphByOrbits ( <G>, <orbits> )
##
##  <Description>
##
##     Returns a connection matrix (on integers) for the coil graph <A>G</A> using
##     the orbits <A>orbits</A>.
##
##  </Description>
##
ConnectionMatrixOfCoilGraphByOrbits := function(G,orbits)
  local m, Zn, totalOrbits, row, col, x, adjacencies, orbit;

  if Union(orbits) <> Vertices(G) then
    Error("orbits is not a partitio of the vertices of G");
  fi;

  totalOrbits := Size(orbits);
  m := EmptyConnectionMatrix(totalOrbits);
  Zn:= ZmodnZ(Size(orbits[1]));

  for row in [1..totalOrbits] do
    for col in [row .. totalOrbits] do
      x := orbits[row][1];
      orbit := orbits[col];

      adjacencies := Adjacency(G,x);

      m[row][col]:= List(Intersection(adjacencies, orbit),i->Position(orbit,i))-1;
      Sort(m[row][col]);

      m[col][row]:= List(-1*(m[row][col]*One(Zn)), Int);
      Sort(m[col][row]);
    od;
  od;

  return m;
end;


CoilCliqueGraph :=function(coilG)

    local kG, gamma, orbits, connectionMatrix,r, cG;

    kG := CliqueGraph(coilG.G);
    gamma := Permutation(coilG.gamma, Cliques(coilG.G), OnSets);
    orbits := Orbits(Group(gamma), Vertices(kG));
    cG := CoilGraph(rec(
      G:=kG,
      verticesPerOrbit:=Size(orbits[1]),
      totalOrbits := Size(orbits),
      orbits := orbits,
      gamma := gamma,
      connectionMatrix := ConnectionMatrixOfCoilGraphByOrbits(kG, orbits)
    ));

    return cG;
end;

InstallOtherMethod(CliqueGraph,[IsRecord],function(coilG)
   ###### INVESTIGAR PORQUE ESTE METODO EN YAGS REGRESA SIEMPRE OBJETOS NO MUTABLES
   return [];
end);


######################################################################
##
#F  CliqueOrbits ( <coil> )
##
##  <Description>
##
##     Given a coil graph, it returns the cliques that define the orbits
##     of its clique graph (wich is a coil graph).Each clique is
##     returned as an array formed by the intersection of the clique with
##     the orbits of the coil graph <A>coil</A>.
##
##  </Description>
##
CliqueOrbits := function(coil)
   local cOrbs, Kcoil, clique, cliqueOrbit, i, j, offset;

   cOrbs := [];
   Kcoil := CoilCliqueGraph(coil);
   for i in [1..Kcoil.totalOrbits] do
      clique := VertexNames(Kcoil.G)[Kcoil.orbits[i][1]];
      cliqueOrbit := [];
      for j in [1..coil.totalOrbits] do
         offset := -1*(j-1)*coil.verticesPerOrbit;
         Add(cliqueOrbit, Intersection(clique, coil.orbits[j])+offset);
      od;
      Add(cOrbs, cliqueOrbit);
   od;

   return cOrbs;
end;


######################################################################
##
#F  AddOffsetToCoilClique ( <n>, <offset>, <q> )
##
##  <Description>
##
##     Returns the coil clique $(q-1)*Zn(n)+offset$. <A>n</A> must be
##     the number of vertices of the coil graph that has the clique
##     coil <q>, which must be a matrix (each row is an orbit of the
##     coil graph).
##
##     Important: the returned clique is on [0..(n-1)].
##
##  </Description>
##
AddOffsetToCoilClique := function(n,offset, q)
  local Zn, qo, orbit;

  Zn := ZmodnZ(n);
  qo := (q-1)*One(Zn) + offset*One(Zn);

  for orbit in [1..Size(q)] do
     qo[orbit] := List(qo[orbit], Int);
     Sort(qo[orbit]);
  od;

  return qo;
end;


######################################################################
##
#F  FindCoilCliqueBoundOffset ( <n>, <maxJump>, <q> )
##
##  <Description>
##
##     Given a coil graph of <A>n</A> vertices per orbit, that has
##     a maximum jump <A>maxJump</A> and a clique <q>, wich must be a
##     matrix (each row contains the vertices of the clique in that orbit),
##     returns a $offset$, if exists, such that:
##
##      $(q-1)*Zn(n)+offset \subset {0..totalOrbits-1}x{0...maxJump}$,
##
##     otherwise, throws an error.
##
##     Important: the clique <q> must be on [1..n] but the offset is
##     computed for q-1 (i.e. the q-1 is a clique on [0..n-1]).
##
##  </Description>
##
FindCoilCliqueBoundOffset := function(n, maxJump, q)
   local Zn, container, boundOffset, offset, isBound, orbitNumber, qi, zqi;

   Zn := ZmodnZ(n);
   container := [0.. maxJump]*One(Zn);

   boundOffset := -1;
   for offset in [0..(n-1)]*One(Zn) do

      isBound := true;
      for orbitNumber in [1..Size(q)] do
         qi:=q[orbitNumber]-1;
         zqi:=qi*One(Zn);
         if IsSubset(container, zqi + offset)=false then
            isBound := false;
            break;
         fi;
      od;

      if isBound=true then
        boundOffset := offset;
        break;
      fi;
   od;

   if boundOffset = -1 then
     Error("Can't center coil clique orbit:", q, ", using n:",n," and max jump:", maxJump, "\n");
   fi;

   return Int(boundOffset);
end;


######################################################################
##
#F  FindCirculantCliqueBoundOffset ( <n>, <maxJump>, <q> )
##
##  <Description>
##
##     Given a circulant graph of <A>n</A> vertices that has
##     a maximum jump <A>maxJump</A> and a clique <q>, wich must be a
##     list, returns an $offset$ in the set [0... n-1] (if exists), such that:
##
##      $(q-1)*Zn(n)+offset \subset \{0...maxJump\}$,
##
##     otherwise throws an error.
##
##     Important: the clique <q> must be on [1..n] but the offset is
##     computed for q-1 (i.e. the q-1 is a clique on [0..n-1]).
##
##  </Description>
##
FindCirculantCliqueBoundOffset := function(n, maxJump, q)
   return FindCoilCliqueBoundOffset(n, maxJump, [q]);
end;



######################################################################
##
#F  CenterCliqueOrbit ( <n>, <maxJump>, <q> )
##
##  <Description>
##
##     If the clique orbit <A>q</A> is bounded
##     by the rectangle R_0(0,maxJump) after aplying some offset, then
##     returns the centered clique orbit:
##     i.e., $q-((max(q)+min(q))/2)$.
##
##     Important: the returned clique is on [0..(n-1)]
##
##  </Description>
##
CenterCliqueOrbit := function(verticesPerOrbit, maxJump, q)
  local offset, qo, qmax, qmin, qi, qimax, qimin, center;

          offset := FindCoilCliqueBoundOffset(verticesPerOrbit, maxJump,q);
          qo := AddOffsetToCoilClique(verticesPerOrbit, offset, q);

          #Print("q-1:", q-1, " has offset:", offset, "\n");
          #Print("(q-1)+offset:", qo, "\n");

          qmax :=0;
          qmin :=verticesPerOrbit;

          for qi in qo do
             if Size(qi)>0 then
                 qimax := Maximum(qi);
                 qimin := Minimum(qi);

                 if qimax > qmax then
                    qmax := qimax;
                 fi;

                 if qimin < qmin then
                    qmin := qimin;
                 fi;

                #Print("qi:",qi,", qimax:", qimax, ", qimin:", qimin, ", qmax:", qmax, ", qmin:", qmin, "\n");
            fi;

          od;

          center := Int((qmax + qmin)/2);
          #Print("(q-1)+offset:",qo, " has center:", center, "\n\n");

          return qo-center;
end;

CenterCliqueOrbit2 := function(n, maxjmp, q)
 local Zn, i, qz, subSet, cOrbs;

 subSet := Concatenation([n-(maxjmp-1)..n], [1..maxjmp+1]);
 Zn := ZmodnZ(n);
 cOrbs := [];
 for i in [0.. n-1] do
    qz := List((q-1)*One(Zn)+i,e->List(e,Int))+1;
    Print("subset:", subSet, "\n");
    Print("qz:",qz, "\n\n");
    if ForAll(qz, sq->IsSubset(subSet,sq)) then
       Add(cOrbs, qz);
       #break;
    fi;
 od;

 return cOrbs;

end;

######################################################################
##
#F  CenteredCliqueOrbits ( <coil>, <maxJump> )
##
##  <Description>
##
##     If every clique orbit of the coil graph <A>coil</A> is bounded
##     by the rectangle R_0(0,maxJump) after aplying some offset, then
##     returns the centered clique orbits:
##     i.e., each clique orbit q is transformed to $q-(max(q)+min(q))/2$.
##
##     Important: the returned centered cliques are on [0..(n-1)]
##
##  </Description>
##
CenteredCliqueOrbits := function(coil, maxJump)
   local cOrbs, centeredcOrbs, i, q, qo, qc, offset, qmax, qmin, qi, qimax, qimin, center, Zn;

   cOrbs := CliqueOrbits(coil);
   Zn := ZmodnZ(coil.verticesPerOrbit);
   centeredcOrbs := [];

   i:=0;
   for q in cOrbs do
      i:=i+1;
      #Print("orbit number:",i,"\n");
      #Print(q,"\n");

      if Size(q) > 0 then
          qc:=CenterCliqueOrbit(coil.verticesPerOrbit, maxJump, q);
          #Print("centered:",qc,"\n");
          Add(centeredcOrbs, qc);
      else
        Add(centeredcOrbs, []);
      fi;
   od;

   return centeredcOrbs;
end;


############################################################################
##
#F  CirculantCliqueOrbitAdjacencies( <qa>, <qb> )
##
##  <Description>
##
##     Returns the adjacencies (a list) of the clique orbit <A>qa</A> with
##     the clique orbit <A>qb</A> for a circulant graph (i.e. T^{ab}, where
##     T is the connection matrix of the clique graph of a circulant graph).
##
##  </Description>
##
CirculantCliqueOrbitAdjacencies := function(qa,qb)
    local Tab, x, y;

    Tab := [];
    for x in qa do
       for y in qb do
          Add(Tab, x-y);
       od;
    od;

    return Union([],Tab); #cheap trick to order and filter repeated values
end;


############################################################################
##
#F  CoilCliqueOrbitAdjacencies( <qa>, <qb> )
##
##  <Description>
##
##     Returns the adjacencies (a matrix) of the clique orbit <A>qa</A> with
##     the clique orbit <A>qb</A> for a coil graph (i.e. T^{ab}, where
##     T is the connection matrix of the clique graph of a coil graph).
##
##     The clique orbits <A>qa</A> and <A>qb</A> must be matrices (with the
##     same number of rows) such that the $i$-th row is a list of the
##     vertices of the $i$-th orbit of the coil graph that belongs to clique
##     orbit (<A>qa</A> or <A>qb</A>).
##
##  </Description>
##
CoilCliqueOrbitAdjacencies := function(qa,qb)
   local Tab, row, circulantClique;

   Tab := [];
   if Size(qa)<>Size(qb) then
      Error(qa, " and ", qb, " aren't the same size");
   fi;

   for row in [1..Size(qa)] do
      Tab := Union(Tab, CirculantCliqueOrbitAdjacencies(qa[row],qb[row]));
   od;

   return Tab;
end;



############################################################################
##
#F  ConnectionMatrixOfCoilGraphByCliqueOrbits( cOrbits )
##
##  <Description>
##
##     Returns the connection matrix of the clique graph of a coil graph
##     with cliques <A>cOrbits</A>: this cliques are supposed to be the
##     clique orbits of the clique graph of a coil graph.
##  </Description>
##
ConnectionMatrixOfCoilGraphByCliqueOrbits := function(cOrbits)
   local m, totalOrbits, row, col;

   totalOrbits := Size(cOrbits);
   m:= EmptyConnectionMatrix(totalOrbits);
   for row in [1..totalOrbits] do
     for col in [row .. totalOrbits] do
        m[row][col] := CoilCliqueOrbitAdjacencies(cOrbits[row],cOrbits[col]);
        m[col][row] := -m[row][col];
        Sort(m[col][row]);
     od;
   od;

   return m;
end;


ConnectionMatrix := function(C)
  return List([1..C.totalOrbits], i-> List([1..C.totalOrbits], j-> Connections2(C,i,j)));
end;

OrbitGraph := function(m, Tad)
  return GraphByRelation([1..Size(m)], function(x,y)
        return m[x][y] in Tad;
      end:
      GraphCategory:=LooplessGraphs
    );
end;

PrintConnectionMatrix := function(C)
  local i,j;

  for i in [1.. C.totalOrbits] do
    Print("{");
    for j in [i..C.totalOrbits] do
      Print(Connections2(C,i,j), "  ");
    od;
    Print("}\n");
  od;
end;



RemoveOrbit := function(m,i)
   local rows, cols;

   rows := Concatenation([1..i-1],[i+1..Size(m)]);
   cols := rows;

   return List(rows, i-> List(cols, j->m[i][j]));
end;


#################################################################################
##                  FOR COIL Cn_13678
#################################################################################

FixCenteredCliqueCoil := function(q)
   local qi, offset;

   offset := 0;
   for qi in q do
      if -2 in qi then
         offset := -1;
         break;
      fi;
   od;

   return q+offset;
end;

FixCenteredCliqueCoils := function(coilsCliques)
   return List(coilsCliques, q->FixCenteredCliqueCoil(q));
end;

IsCliqueComposedOfStandardSmallCliques := function(q)
   local isOk, catalog, qi, i;

   isOk:= true;
   catalog := [
     [], [0], [-4], [4], [-3,3], [-3,0,3], [-3,3,4], [-4,-3,3], [-4,-3,3,4]
   ];

   for i in [1..Size(q)] do
     qi := q[i];
     if (qi in catalog)=false then
        Print("clique part with index:", i, " of clique orbit is not standard:\n");
        Print(qi, "\n");
        isOk := false;
        break;
     fi;
   od;

   return isOk;
end;

IsListOfCliquesComposedOfStandardSmallCliques := function(list)
  local isOk, q, i;

  isOk := true;
  for i in [1..Size(list)] do

     q := list[i];
     if IsCliqueComposedOfStandardSmallCliques(q)=false then
        Print("clique orbit with index:", i, " is not standard: \n");
        Print(q, "\n");
        isOk := false;
        break;
     fi;
  od;

  return isOk;
end;

IsCliqueStandardOrbit := function(clique)
    local adj, stdOrbits;

    stdOrbits := List([[0,1,3,6,7,8],[0,1,6,7,8], [0,1,6,7], [0,6], [0,3,6]], r-> Set(Concatenation(r,-r)));
    adj := CoilCliqueOrbitAdjacencies(clique,clique);

    if adj in stdOrbits then
      return true;
    else
      return false;
    fi;
end;

IsListOfCliquesStandardOrbit := function(list)
   local isStd, i, clique;

   isStd := true;

   for i in [1.. Size(list)] do
      clique := list[i];
      if not IsCliqueStandardOrbit(clique) then

        Print("clique ", i, " is not standard\n");
        isStd := false;
        break;
      fi;
   od;

   return isStd;
end;



##############test#############
#K:=CliqueGraph;
#n:=13;
#C:=Circulant(n,[1,2,4]);; KC:=CliqueGraph(C);; K2C:=CliqueGraph(KC);
#Zn:=ZmodnZ(n);
#m:=[ [ [1,2,4], [-3,-2,-1,0,1,2,3]], [ [-3,-2,-1,0,1,2,3],[1,2]] ];;
#r:=CoilGraphByConnectionMatrix(Zn, m);
#Print(IsIsomorphicGraph(r.G,K2C), "\n");

#PrintListOneItemPerLine(ConnectionMatrixOfCoilGraphByOrbits(r.G, r.orbits));

#co:=CoilGraphByConnectionMatrix(Zn,[[[1,2,4]] ]);
#Print(IsIsomorphicGraph(co.G,C), "\n");
#Kco := CoilCliqueGraph(co);
#Print(IsIsomorphicGraph(CoilGraphByConnectionMatrix(Zn, Kco.connectionMatrix).G,KC), "\n");

#K2co := CoilCliqueGraph(Kco);
#Print(IsIsomorphicGraph(CoilGraphByConnectionMatrix(Zn, K2co.connectionMatrix).G,K2C), "\n");

#n:=15;
#Zn:=ZmodnZ(n);
#co:=CoilGraphByConnectionMatrix(Zn,[[[1,2,4]] ]);

#
#IMPORTANTE: GAP es tan estupido que si uno imprime o guarda texto en un archivo (incluso usando output stream)
#            escapa el caracter especial \ para indicar que la linea continua, esto es un problema al compilar
#            archivos de GAP
#            ¿¡¡QUE DIABLOS PENSABAN LOS CREADORES DE GAP!!? ¿¿¿ACASO SON ESTUPIDOS????
#            por eso tengo que poner \n
#
ListToTokenizedString := function(delimiters, list)
   local str, item, delimiter, counter, newd;

  #@TODO: check that the dimension of the list is equals to the size of delimiters list
  str := "";

  if IsList(list) then
      delimiter := delimiters[1];
      if Size(list) = 0 then
        str := "\\EMPTY_SET_CODE\n";
      else
        counter := 1;
        for item in list do
          #Print("list:", list, "\n");
          #Print("del:",delimiters,"\n");
          newd := delimiters{[2..Size(delimiters)]};
          #Print("newd:",newd,"\n");
          str := Concatenation(str, ListToTokenizedString(newd,item));
          if counter < Size(list) then
             str := Concatenation(str, delimiter);
          fi;
          counter := counter + 1;
        od;
      fi;
  else
      str := Concatenation(String(list),"\n");
  fi;

  return str;
end;
#gap> ListToTokenizedString(["\\\\",";",","], [[[-4,-3,3],[-3,0,3],[0,1]],[[-3,3],[],[-4,3,3,4]],[[0],[0],[-4]]]);
#"-4,-3,3;-3,0,3;0,1\\\\-3,3;\\EMPTY_SET_CODE;-4,3,3,4\\\\0;0;-4"
# ListToTokenizedString(["\\\\",";",","], q);

CliqueOrbitsToTikz := function(cliqueOrbits)
   local str,tikzTempl;
   str := ListToTokenizedString(["\\\\",";",","], cliqueOrbits);
   tikzTempl := ReplacedString(ReadAll(InputTextFile("latex/tikz/clique_orbits2.tmpl")), "TOKEN_BIG_CLIQUES", str);
   return tikzTempl;
end;
# CliqueOrbitsToTikz(q);

CliqueOrbitsToLatex := function(title,cliqueOrbits)
   local str;

   str := CliqueOrbitsToTikz(cliqueOrbits);
   str := ReplacedString(ReadAll(InputTextFile("latex/tikz/pdf-template.tmpl")), "<>", str);
   str := ReplacedString(str, "TITLE_PDF", title);

   return str;
end;
# CliqueOrbitsToLatex("Clanes de $K^{1}(C_{25}(1,3,6,7,8))$", q);

CliqueOrbitsToPDF := function(dir, filename, title, cliqueOrbits)
   local pdfTempl, cmd, filepath, output;

   pdfTempl := CliqueOrbitsToLatex(title, cliqueOrbits);
   filepath := Concatenation(dir,"/",filename);
   #output := OutputTextFile( filepath, true );;
   #AppendTo( output, pdfTempl);
   #CloseStream(output);
   PrintTo(filepath, pdfTempl);
   cmd := Concatenation("lualatex --interaction=nonstopmode -output-directory ",dir, " ", filepath);
   cmd := Concatenation(cmd," | awk 'BEGIN{IGNORECASE = 1}/warning|!/,/^$/;'");
   #Print(cmd);
   Exec(cmd);;
end;

#CliqueOrbitsToPDF("latex","test.tex","Clanes de $K^{1}(C_{25}(1,3,6,7,8))$", q);
#CliqueOrbitsToPDF("latex","coils.tex","Clanes de $K^{1}(C_{25}(1,3,6,7,8))$", q);


CanonicalSmallCliques := function()
   return [
      #[],
      [-4],
      [4],
      [-3,3],
      [-3,0,3],
      [-3,3,4],
      [-4,-3,3],
      [-4,-3,3,4]
  ];
end;

CanonicalJumps := function()
  local jumps1, jumps2;

  jumps1 := [
     [1,3,6,7,8], #domos y carpas
     [3,6],       #carpas
     [1,6,7,8],   #domos
     [1,6,7],     #gotas
     [6]          #arcos
  ];

  #orden en libreta
  jumps2 := [
     [1,3,6,7,8],
     [1,6,7,8],
     [1,6,7],
     [6],
     [3,6]
  ];

  return jumps1;
end;

CanonicalJumps2 := function()

   return List(CanonicalJumps(), J-> Set(Concatenation(J,-1*J,[0])));
end;

CanonicalJumps3 := function(n)
   local canOrbits, canOrb;

   canOrbits := CanonicalJumps2()*One(ZmodnZ(n));
   canOrbits := List(canOrbits, canOrb->List(canOrb, i-> Int(i)));

   for canOrb in canOrbits do
     Sort(canOrb);
   od;

   return canOrbits;
end;


OrbitCount := function(coils)
   local n, coilsStat, table, canJumps, canJumpIndex, canJump, jumpStat, coilStat, pos;

   n := coils[1].verticesPerOrbit;
   coilsStat :=  List(coils, coil ->
               Collected( List([1..Size(coil.connectionMatrix)], i->coil.connectionMatrix[i][i]))
          );

   table := [];
   canJumps := CanonicalJumps3(n);

   for canJumpIndex in [1 .. Size(canJumps)] do
      canJump := canJumps[canJumpIndex];
      jumpStat := [];

      for coilStat in coilsStat do
          pos := Position( List(coilStat, i-> i[1]), canJump);
          if pos=fail then
            Add(jumpStat, 0);
          else
            Add(jumpStat, coilStat[pos][2]);
          fi;
      od;

      table[canJumpIndex] := jumpStat;
   od;

   return table;
end;


OrbitCountTableToPDF := function()
   local rowHeaders, colHeaders;

   rowHeaders := List([1..Size(CanonicalJumps())],j-> Concatenation("$C_", String(j-1), "$"));
   Print(rowHeaders);
   colHeaders := CanonicalJumps();

end;


SmallCliquesToTablePDF := function()
  local str, cliStr, i,j, smallCliques, smallClique, adj, cliqueA, cliqueB, rowSpace;

  smallCliques := CanonicalSmallCliques();

  rowSpace := "";
  for i in [1..Size(smallCliques)] do
     rowSpace := Concatenation(rowSpace, "&");
  od;
  rowSpace := Concatenation(rowSpace, "\\\\\n");


  str := String("");
  for i in [1..Size(smallCliques)] do
     str := Concatenation(str, "|c");
  od;
  str := Concatenation(str, "|");

  str := Concatenation("\n\n\\begin{tabular}{|c", str ," } \n");
  str := Concatenation(str, "\\hline\n");
  str := Concatenation(str, "&\n");



  ## headers
  for i in [1..Size(smallCliques)] do

    smallClique := smallCliques[i];
    str := Concatenation(str, "\\begin{tikzpicture}\n");
    cliStr := ReplacedString(String(smallClique),"[","{");
    cliStr := ReplacedString(cliStr,"]","}");
    str := Concatenation(str, "\\draw pic {orbAndSmallClique=0/0/0.5/", cliStr, "/", String(Size(smallClique)), "};\n");
    str := Concatenation(str, "\\end{tikzpicture}\n");

    if i = Size(smallCliques) then
      str := Concatenation(str, "\\\\\n");
    else
      str := Concatenation(str, "&\n");
    fi;
  od;
  str := Concatenation(str, "\\hline\n");

  ## body
  for i in [1..Size(smallCliques)] do

    str := Concatenation(str, rowSpace);
    cliqueA := smallCliques[i];
    str := Concatenation(str, "\\begin{tikzpicture}\n");
    cliStr := ReplacedString(String(cliqueA),"[","{");
    cliStr := ReplacedString(cliStr,"]","}");
    str := Concatenation(str, "\\draw pic {orbAndSmallClique=0/0/0.5/", cliStr, "/", String(Size(cliqueA)), "};\n");
    str := Concatenation(str, "\\end{tikzpicture}\n");
    str := Concatenation(str, "&\n");

    for j in [1..Size(smallCliques)] do
      cliqueB := smallCliques[j];
      adj :=  CirculantCliqueOrbitAdjacencies(cliqueA,cliqueB);
      str := Concatenation(str, " ", String(adj), " \n");

      if j = Size(smallCliques) then
        str := Concatenation(str, "\\\\\n");
      else
        str := Concatenation(str, "&\n");
      fi;
    od;
    str := Concatenation(str, rowSpace);
    str := Concatenation(str, rowSpace);
    str := Concatenation(str, "\\hline\n");

  od;

  str := Concatenation(str, "\\hline\n");
  str := Concatenation(str, "\\end{tabular}\n");
  str := Concatenation(str, "\n\n");

  str := ReplacedString(ReadAll(InputTextFile("latex/tikz/pdf-template-table.tmpl")), "<>", str);
  PrintTo("latex/tables/tabla-clancitos.tex", str);

end;


ClanesPiezaAPDF := function(n, piezaA, piezaB)
  local nombreA, nombreB, nombreArch, clOrbits, coil, q, clanA, clanB;

  nombreA := piezaA[1];
  clanA := piezaA[2];

  nombreB := piezaB[1];
  clanB := piezaB[2];

  Print("construyendo bobina (", nombreA , ", ", nombreB, ")\n");
  nombreArch := Concatenation(nombreA,"_", nombreB,".tex");
  clOrbits := [[clanA], [clanB]];
  #Print("clanes orbita:\n");
  #PrintListOneItemPerLine(clOrbits);
  #Print("\n");
  coil := CoilGraphByConnectionMatrix(n, ConnectionMatrixOfCoilGraphByCliqueOrbits(clOrbits));
  q:=CenteredCliqueOrbits(coil,8);;
  q:=FixCenteredCliqueCoils(q);;
  #Print("clanes centrados:\n");
  #PrintListOneItemPerLine(q);
  #Print("\n");

  #if IsListOfCliquesComposedOfStandardSmallCliques(q)=false or IsListOfCliquesStandardOrbit(q)=false then
  #  Error("pieza:", piezaA[1], " y pieza:", piezaB[1], " producen clanes orbita que no son estandard!\n");
  #fi;

  CliqueOrbitsToPDF("latex/piezas", nombreArch, " piezas ", q);
end;


ReportarPiezasAPDF := function()
  local carpa, domo, gotIz, gotDer, arco, pIzq, pDer, pCen, piezas, piezaA, piezaB, n;

  n:=25;
  carpa:=[-3,0,3];
  domo := [-4,-3,3,4];
  gotIz := [-4,-3,3];
  gotDer := [-3,3,4];
  arco := [-3,3];
  pIzq := [-4];
  pDer := [4];
  pCen := [0];

  piezas := [
    ["domo", domo],
    ["carpa",carpa],
    ["gotIz", gotIz],
    ["gotDer", gotDer],
    ["arco", arco],
    ["pIzq", pIzq],
    ["pDer", pDer],
    ["pCen", pCen]
  ];

  for piezaA in piezas do
     for piezaB in piezas do
        ClanesPiezaAPDF(n, piezaA, piezaB);
     od;
  od;

end;


ReportarTablaClancitos := function()
    local clancitos;

    clancitos := [
      #[],
      [-4],
      [4],
      [-3,3],
      [-3,0,3],
      [-3,3,4],
      [-4,-3,3],
      [-4,-3,3,4]
    ];

    PrintListOneItemPerLine(List(Cartesian(clancitos,clancitos), r->CoilCliqueOrbitAdjacencies([r[1]],[r[2]])));

end;


Report2CoilToPdf := function (q,k,i,j)
  local qa,qb,T,qt,fiNam;

  qa := q[i];
  qb := q[j];
  T:=CoilGraphByConnectionMatrix(25, ConnectionMatrixOfCoilGraphByCliqueOrbits([qa, qb]));;
  qt := CenteredCliqueOrbits(T,8);;
  qt := FixCenteredCliqueCoils(qt);;
  fiNam := Concatenation("clanes_q", String(k), "[", String(i), ",", String(j), "].tex");
  Print("generating:", fiNam, " \n");
  CliqueOrbitsToPDF("latex",fiNam," ", qt);
end;


SortCenteredCliqueOrbits := function(list)
  local jumpsOrder, origJumps, i, jump, pos, jumpsPos;

  jumpsOrder := CanonicalJumps();

  jumpsOrder := List(jumpsOrder, J-> Set(Concatenation(J,-1*J,[0])));
  origJumps := List(list,q-> CoilCliqueOrbitAdjacencies(q,q));

  #Print("original jumps: \n");
  #PrintListOneItemPerLine(origJumps);

  jumpsPos := [];
  for i in [1..Size(origJumps)] do
     jump := origJumps[i];
     pos := Position(jumpsOrder, jump);
     Add(jumpsPos, pos);

     if pos = fail then
        Error("clique orbit ", i, " produces an unrecognized jump:", jump, " \n");
     fi;
  od;

  SortParallel(jumpsPos, list);
end;
