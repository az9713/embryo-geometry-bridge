"""Separate high-precision replay of saved angle signs and slab bound constants."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json
from flint import arb,ctx

HERE=Path(__file__).resolve().parent


def main():
    ctx.prec=320
    p=HERE/'slab_verification.json';source=json.loads(p.read_text())
    assert source['written_theorem_sha256']==sha256((HERE/'THEOREM.md').read_bytes()).hexdigest()
    A=lambda q:arb(str(q))
    angles=[]
    for row in source['interval_angle_examples']:
        q=Q(row['compression_ratio_exact']);d=q/4;lo,hi=map(Q,row['angle_bracket_exact'])
        mean=lambda t:(1-A(d))*A(t).cos()+A(d)*A(t).sin()/A(t)
        assert (mean(lo)-A(q)).lower()>0
        assert (mean(hi)-A(q)).upper()<0
        assert (A(hi).cos()-A(q/2)).lower()>0
        C=Q(row['collar_C_exact']);assert C>=18/q
        assert q/2-Q(9,2)/C>=q/4
        assert Q(row['body_J_lower_exact'])==q/4
        angles.append({'compression_ratio':str(q),'both_mean_signs_and_tangent_guard':True})
    q=Q(4,5);C=18/q;M=Q(3,2)
    row=next(x for x in source['interval_angle_examples'] if Q(x['compression_ratio_exact'])==q)
    hi=Q(row['angle_bracket_exact'][1])
    # Rebuild the compact set, not the saved decimal constants.
    norm2=Q(3)+(2*M/C)**2+(1+M/C)**2
    Jlo=q/4;Jhi=1+q/4
    logs=[abs(A(Jlo).log()).upper(),abs(A(Jhi).log()).upper()]
    L=max(logs)
    MW=(A(norm2)+3+2*L)/2+L*L/2+A(norm2)/2
    energy_coefficient=2*A(C)*MW+8*A(hi)**2/A(q/4)
    rows=[]
    for r in source['illustrative_bound_rows']:
        e=Q(1,r['N']);assert Q(r['epsilon_exact'])==e and Q(r['b_exact'])==e**3
        assert Q(r['collar_depth_exact'])==C*e and 2*C*e<1
        # b=e^3 makes the bound exactly coefficient*e for this sequence.
        energy=energy_coefficient*A(e);position=A(5).sqrt()*A(e)
        assert Q(str(energy.upper().fmpq()))<=Q(r['energy_error_upper_exact'])
        assert Q(str(position.upper().fmpq()))<=Q(r['uniform_position_error_upper_exact'])
        rows.append({'N':r['N'],'energy_and_position_upper_bounds_confirmed':True})
    out={'status':'PASS_SEPARATE_SLAB_CONSTANT_REPLAY','source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'producer_receipt_sha256':sha256(p.read_bytes()).hexdigest(),
         'precision_bits':320,'angle_checks':angles,'bound_checks':rows,
         'scope':'Separate constant/sign arithmetic; written analytic proof still requires mathematical assessment. Not an external peer review.'}
    (HERE/'slab_independent_constants.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':out['status'],'angle_checks':len(angles),'bound_rows':len(rows)}))


if __name__=='__main__':main()
