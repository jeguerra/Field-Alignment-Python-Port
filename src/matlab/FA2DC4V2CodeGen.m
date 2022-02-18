a0 = coder.typeof(0,[inf,inf,inf]);
a1 = coder.typeof(0,[inf,inf]);
a2 = coder.typeof(0,[1,1]);
codegen -config:lib FA2DC4V2LIB.m -args {a0,a0,a0,a0,a0,a1,a1,a1,a2,a2,a2,a2,a2,a2};
codegen -config:dll FA2DC4V2LIB.m -args {a0,a0,a0,a0,a0,a1,a1,a1,a2,a2,a2,a2,a2,a2};
codegen -config:mex FA2DC4V2LIB.m -args {a0,a0,a0,a0,a0,a1,a1,a1,a2,a2,a2,a2,a2,a2};

