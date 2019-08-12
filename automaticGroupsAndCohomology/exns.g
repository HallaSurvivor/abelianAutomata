# GAP files for working with automatic group extensions
#
# Chris Grossack, 2019

# Syntax Reminders:
# SymmetricGroup(n);
# DihedralGroup(n); (group of size n)


######################################################
# CheckExn(K,G,Q)                                    #
#                                                    #
# true iff G has K a normal subgroup with quotient Q #
######################################################

CheckExn := function(K,G,Q)
  local ks, phi, H;

  # all (conjugacy classes of) injections K to G
  ks := IsomorphicSubgroups(G,K); 
  for phi in ks do
    H := Image(phi);
    if IsNormal(G,H) then
      if IdGroup(FactorGroup(G,H)) = IdGroup(Q) then
        return true;
      fi;
    fi;
  od;
  return false;
end;


##########################################################
# DoesExnExist(K,Q)                                      #
#                                                        #
# return a list of all G such that 1 -> K -> G -> Q -> 1 #
##########################################################

DoesExnExist := function(K,Q)
  local g, n, ret, i, G;

  g := Size(K) * Size(Q);
  n := NumberSmallGroups(g);
  ret := [];
  for i in [1 .. n] do
    Print(i); Print("/"); Print(n); Print("\n");
    G := SmallGroup(g, i);
    if CheckExn(K,G,Q) then Add(ret, (g,i)); fi;
  od;
  return ret;
end;


#############################################################
# Simplify(G)                                               #
#                                                           #
# Print the generators/relations of a group isomorphic to G #
# Return the new group                                      #
#############################################################

Simplify := function(G)
  local G2;
  G2 := SimplifiedFpGroup(Image(IsomorphismFpGroup(G)));
  Print(G2);
  Print("\n");
  Print(RelatorsOfFpGroup(G2));
  Print("\n");
  return G2;
end;
