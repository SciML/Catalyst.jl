# Re-exported upstream API
```@meta
CurrentModule = Catalyst
```

Catalyst re-exports selected public APIs from upstream SciML and symbolic
computing packages for modeling convenience. The source-site docstrings for
these names are rendered below from their defining modules.

## ModelingToolkitBase
```@autodocs; canonical=false
Modules = [ModelingToolkitBase]
Public = true
Private = false
```

## Symbolics
```@autodocs; canonical=false
Modules = [Symbolics]
Public = true
Private = false
```

## SymbolicUtils
```@autodocs; canonical=false
Modules = [SymbolicUtils, SymbolicUtils.Code]
Public = true
Private = false
```

## SciMLBase
```@autodocs; canonical=false
Modules = [SciMLBase]
Public = true
Private = false
```

## JumpProcesses
```@autodocs; canonical=false
Modules = [JumpProcesses]
Public = true
Private = false
```

## Other re-exports
```@meta
CurrentModule = Main
```

```@docs; canonical=false
BipartiteGraphs.BipartiteGraph
CommonSolve.solve
TermInterface.arguments
TermInterface.iscall
TermInterface.operation
TermInterface.sorted_arguments
UnPack.@pack!
UnPack.@unpack
```

```@autodocs; canonical=false
Modules = [UnPack]
Public = true
Private = false
```
