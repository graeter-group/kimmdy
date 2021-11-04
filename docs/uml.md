# UML

(Re)render the diagrams with

```bash
plantuml -tsvg uml.md
```

or preview online with [plantuml](https://plantuml.com).

## Classes

```
@startuml classuml
class           RunManager{
config : Config
config_reactions
state : Enum
MDparams
rates 
run_md_min()
run_md_eq()
run_md()
query_reactions()
make_decision()
run_recipe(ConversionRecipe)
next()
}
class          Config{
path
...
read_yaml()
}
class           MDManager{
write_mdp()
do_md()
}
class           ConversionRecipe{
type
atom_idx
}
class           DecisionStrategy{
}
class           Reaction{
recipe: ConversionRecipe
get_rates()
}
class           HAT{

}
class           CoordinateChanger{
}
class           TopologyChanger{
}
package GROMACS{
interface GROMACSCli{
}
}

Reaction --|> HAT
Config -- RunManager : reads <
RunManager - MDManager : > starts
MDManager - GROMACSCli : > interacts
RunManager -- Reaction : > starts
RunManager -- DecisionStrategy : > starts
ConversionRecipe - Reaction: < generates
CoordinateChanger - RunManager: < starts
TopologyChanger - RunManager: < starts
TopologyChanger -[hidden]-> CoordinateChanger
MDManager -[hidden]-> DecisionStrategy
@enduml
```

![](classuml.svg)

```
@startuml flow
start
:Read .yaml;
:Run MD_eq;
repeat
 :Run MD_prod;
 :Query reactions;
 :Make reaction decision;
 :Change system coordinates;
 :Change system topology;
repeatwhile (itermax reached) is (no)
-> yes;
stop
@enduml
```

![](flow.svg)

[link](www.plantuml.com/plantuml/png/bLB1Zfim4BtFL_YuUo0VW5RHhaZLNb0bTczLXImyIwnWczYc3QBzzmuaO0kuDBdC36zctdlZL-UvzPrge6guSopyYaxdNFCQxG2LqP-oPYdBfk2HbnPvvQNH3cYAH_h-HNSAybFBBLSEB1KT0zlfKebIIVtqF2TuNM8AhXtQFeoZYk8NB0LMqaapjrbAMtmY3h_GZlLYAZo3nfidpD-rXZlR0Lhkpt0u780sYBBdgjb3i_oq2FvjfKVYrX9G60fs6zPC1l1zYy2zKQKKjvsEqqFHkn-zgVjX1rCyR1ZWBZZTx84QVJcPkqlhszl70BjqZHLIKrzsvffqxct_CArfJDr7a9PN5xA5VIs-vs_P-m1IUxIVl5fAMIC9I7-OoRCa-NCScS3z29H7BugbR3o5OoyG5PDm0G8SsVGq7OHY4ksR41CHNX4e7fCi5eOnQyJw435oRJSB0rFsmchKE9aF6qDPB9AmyyHWwl_1Vfob4AX_9iVjyM9V0cuZ35vHYqur_m00)

[link](www.plantuml.com/plantuml/png/POz12i8m44NtESKisuKUm88KTDk5Na26PcY3IIPcfcXkRrge6tV_lmzlc5uKiox1cosOqvtGmh5Wy5qjIuJX-g1NPp8bGZMmivJPThNU5ie5Ck6eZgEiQC0d_GXO6ftKi2wN6UD484MK0epCsRg8Il8_AYVsF9Nzydjsdg1nIZdWPzFFP5jm0atarXpEK5QFA2VJKxJrfoy0)
