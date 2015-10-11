# DUNE_ISTL_CG
Conjugate Gradient Methods for DUNE.

Structure is as follows. There is a generic iterative method 

<code> 
template &lt;class Step,class TerminationCriterion&gt; class GenericIterativeMethod;
</code>

which serves as building block for different conjugate gradient methods. Currently there exist two termination criteria

<code>namespace Dune{ namespace KrylovTerminationCriterion{</code>

<code>    template &lt;class real_type&gt; class ResidualBased; </code>
    
<code>    template &lt;class real_type&gt; class RelativeEnergyError;</code>

<code>} }</code>

The step computation is again decomposed into different substeps that work on a common data structure. The general structure is as follows (though most of the steps can be replaced with whatever you want it to be):
 - Application of the preconditioner
 - Computation of the search direction
 - Computation of the scaling parameter for the search direction
 - Treat nonconvexity
 - Update iterate
 - Adjust other data

Based on this GenericStep, different conjugate gradient solvers and the Chebyshev semi-iteration are implemented
The solvers are currently called MyCGSolver, TCGSolver, RCGSolver, TRCGSolver and support different terminatin criteria. 
The syntax is as previously with additional optional template parameter for the termination criterion.
The simplest ways to generate a cg solver(in namespace Dune) are:

<code>auto cg   = make_cg&lt;MyCGSolver,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp);</code>

<code>auto tcg  = make_cg&lt;TCGSolver,KrylovTerminationCriterion::RelativeEnergyError&gt;(A,P,sp);</code>

<code>auto rcg  = make_cg&lt;RCGSolver,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp);</code>

<code>auto trcg = make_cg&lt;TRCGSolver,KrylovTerminationCriterion::RelativeEnergyError&gt;(A,P,sp);</code>

or

<code>auto cg   = MyCGSolver&lt;Domain,Range&gt;(A,P,sp,terminationCriterion);</code>

<code>auto tcg  = TCGSolver&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion);</code>

<code>auto rcg  = RCGSolver&lt;Domain,Range&gt;(A,P,sp,terminationCriterion);</code>

<code>auto trcg = TRCGSolver&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion);</code>
