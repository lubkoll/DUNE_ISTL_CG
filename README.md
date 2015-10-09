# DUNE_ISTL_CG
Conjugate Gradient Methods for DUNE.

Structure is as follows. There is a generic iterative method 

<code>
template &lt;class Step,
          class TerminationCriterion&gt;
class GenericIterativeMethod;
</code>

which serves as building block for different conjugate gradient methods. Currently there exist two termination criteria

<code>
namespace Dune{

namespace KrylovTerminationCriterion{
    template &lt;class real_type&gt;
    class ResidualBased;
    
    template &lt;class real_type&gt;
    class RelativeEnergyError;
  }
}
</code>
Based on a general CGStepImpl a policy-based approach is used to adjust the implementation for truncated, regularized and truncated regularized conjugate gradient methods.
The solvers are called CG, TCG, RCG, TRCG and support different terminatin criteria. 
The simplest way to generate a cg solver is:

<code>
auto cg   = make_cg<CG,KrylovTerminationCriterion::ResidualBased>(A,P,sp);

auto tcg  = make_cg<TCG,KrylovTerminationCriterion::RelativeEnergyError>(A,P,sp);

auto rcg  = make_cg<RCG,KrylovTerminationCriterion::ResidualBased>(A,P,sp);

auto trcg = make_cg<TRCG,KrylovTerminationCriterion::RelativeEnergyError>(A,P,sp);
</code>
