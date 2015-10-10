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

Based on a general CGStepImpl a policy-based approach is used to adjust the implementation for truncated, regularized and truncated regularized conjugate gradient methods.
The solvers are called CG, TCG, RCG, TRCG and support different terminatin criteria. 
The simplest ways to generate a cg solver(in namespace Dune) are:

<code>auto cg   = make_cg&lt;CG,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp);</code>

<code>auto tcg  = make_cg&lt;TCG,KrylovTerminationCriterion::RelativeEnergyError&gt;(A,P,sp);</code>

<code>auto rcg  = make_cg&lt;RCG,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp);</code>

<code>auto trcg = make_cg&lt;TRCG,KrylovTerminationCriterion::RelativeEnergyError&gt;(A,P,sp);</code>

or

<code>auto cg   = CG&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion,maxSteps);</code>

<code>auto tcg  = TCG&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion,maxSteps);</code>

<code>auto rcg  = RCG&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion,maxSteps);</code>

<code>auto trcg = TRCG&lt;Domain,Range,KrylovTerminationCriterion::ResidualBased&gt;(A,P,sp,terminationCriterion,maxSteps);</code>
