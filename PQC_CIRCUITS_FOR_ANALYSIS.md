# Step-5 learned circuits — exact spec for analytic analysis

Convention: 5 qubits, wire 0 = ancilla `a`, wires 1,2 = kept register A (A1,A2), wires 3,4 = discarded register B (B1,B2).  Single-qubit gates are ideal; the noise model puts a 2-qubit depolarizing channel `(1-e2)rho + e2 I/4` after **each CNOT**.  Rotation convention `RX(t)=exp(-i t X/2)` etc.; `Rot(phi,theta,omega) = RZ(omega) RY(theta) RZ(phi)`.  Input `rho_eps = (1-eps)|Phi+><Phi+| + eps I/4` on register (A,B) pairs (A1B1 and A2B2 each carry one noisy Bell copy).

Two representations are given per circuit: (1) the **exact primitive gate list** as compiled (RX/RY/RZ per wire, CNOTs), fully reproducible; (2) the **merged form** — consecutive single-qubit layers composed into one net SU(2) per wire between CNOTs — which is the minimal-depth circuit and the convenient object for a pen-and-paper conjugation of the observables through the CNOTs.  Both implement the identical unitary.


==============================================================================
## 5a / 5b — 14-CNOT full-unitary circuit
==============================================================================

This single circuit compiles the full 5-qubit gadget unitary U = H_a . CSWAP(0;1,3) . CSWAP(0;2,4) . H_a to machine precision, and therefore also realizes the ancilla-|0> isometry exactly (5b). It is the greedy-pruning floor of the gadget-matched ansatz (13 CNOTs unreachable). Under per-CNOT depolarizing noise it defines the operational threshold analysed in pqc_ring_threshold.py.

- CNOT sequence (in order): [(0, 1), (1, 3), (0, 4), (2, 4), (0, 3), (1, 3), (0, 4), (2, 4), (0, 1), (0, 3), (1, 3), (0, 2), (0, 4), (2, 4)]
- CNOT count: 14    rotation blocks: 15
- verification: |Tr(U_target^dag V)|/32 = 1.000000000000  (=1 => exact compilation of U)
- drawn(merged)==solution overlap = 1.000000000000

### (1) Exact primitive gate list
```
    RX(theta=+1.4930918042, wire=0)
    RY(theta=-1.5707963180, wire=0)
    RZ(theta=-1.5390720903, wire=0)
    RX(theta=-0.2969550082, wire=1)
    RY(theta=+0.9823921611, wire=1)
    RZ(theta=+1.7724349797, wire=1)
    RX(theta=+2.6778036789, wire=2)
    RY(theta=+1.7141840304, wire=2)
    RZ(theta=+2.5960971313, wire=2)
    RX(theta=+0.4454536645, wire=3)
    RY(theta=+1.9997639366, wire=3)
    RZ(theta=+0.8912452771, wire=3)
    RX(theta=+1.7669199301, wire=4)
    RY(theta=-1.0504459413, wire=4)
    RZ(theta=+0.5873214048, wire=4)
  CNOT(control=0, target=1)
    RX(theta=+1.5713662174, wire=0)
    RY(theta=-0.6103666788, wire=0)
    RZ(theta=+3.0431100216, wire=0)
    RX(theta=+0.9360362825, wire=1)
    RY(theta=+0.7252233058, wire=1)
    RZ(theta=-2.5099793370, wire=1)
    RX(theta=-0.6253164294, wire=2)
    RY(theta=+1.1932266122, wire=2)
    RZ(theta=+1.4778931020, wire=2)
    RX(theta=-2.0025313187, wire=3)
    RY(theta=-1.5399750521, wire=3)
    RZ(theta=+1.3204525535, wire=3)
    RX(theta=+1.7293807818, wire=4)
    RY(theta=-0.0983673650, wire=4)
    RZ(theta=+1.6633485772, wire=4)
    RX(theta=-2.3968263065, wire=0)
    RY(theta=+2.0700479371, wire=0)
    RZ(theta=-1.5825798616, wire=0)
    RX(theta=+1.1972510484, wire=1)
    RY(theta=-2.3168752670, wire=1)
    RZ(theta=-1.2942813150, wire=1)
    RX(theta=-0.2275212803, wire=2)
    RY(theta=+1.0955969049, wire=2)
    RZ(theta=+2.6381754666, wire=2)
    RX(theta=-1.9247433215, wire=3)
    RY(theta=+0.3310685719, wire=3)
    RZ(theta=-1.2801002672, wire=3)
    RX(theta=+0.2857300053, wire=4)
    RY(theta=+0.2647517124, wire=4)
    RZ(theta=-1.7331655752, wire=4)
  CNOT(control=1, target=3)
    RX(theta=+2.0033977880, wire=0)
    RY(theta=-1.4118223587, wire=0)
    RZ(theta=+3.3368189686, wire=0)
    RX(theta=-2.9009306504, wire=1)
    RY(theta=-0.5434523634, wire=1)
    RZ(theta=-2.3055972826, wire=1)
    RX(theta=-0.7099618076, wire=2)
    RY(theta=-0.9356443357, wire=2)
    RZ(theta=+1.3761679846, wire=2)
    RX(theta=-0.5226719832, wire=3)
    RY(theta=+2.2310891025, wire=3)
    RZ(theta=+2.7300955480, wire=3)
    RX(theta=+0.1652880596, wire=4)
    RY(theta=+2.1258041046, wire=4)
    RZ(theta=-2.9807095096, wire=4)
    RX(theta=-0.7488172002, wire=0)
    RY(theta=+0.5136349785, wire=0)
    RZ(theta=-1.6243643175, wire=0)
    RX(theta=+2.1991227956, wire=1)
    RY(theta=+2.7836384037, wire=1)
    RZ(theta=-2.7104368117, wire=1)
    RX(theta=+2.1646013562, wire=2)
    RY(theta=+0.2387367829, wire=2)
    RZ(theta=+1.0094095319, wire=2)
    RX(theta=-2.9116850730, wire=3)
    RY(theta=+3.0201888903, wire=3)
    RZ(theta=-0.5595264986, wire=3)
    RX(theta=-2.8272167468, wire=4)
    RY(theta=+2.3840547911, wire=4)
    RZ(theta=+1.0638568779, wire=4)
  CNOT(control=0, target=4)
    RX(theta=-1.6041135510, wire=0)
    RY(theta=-1.7901942252, wire=0)
    RZ(theta=+0.9440300208, wire=0)
    RX(theta=-0.6249836150, wire=1)
    RY(theta=-1.4568455372, wire=1)
    RZ(theta=+0.1437169078, wire=1)
    RX(theta=-1.1524037785, wire=2)
    RY(theta=+1.8677356562, wire=2)
    RZ(theta=+1.3557910983, wire=2)
    RX(theta=-1.5487703244, wire=3)
    RY(theta=-0.6808635270, wire=3)
    RZ(theta=-1.4693547253, wire=3)
    RX(theta=-0.6290258796, wire=4)
    RY(theta=+1.5707963420, wire=4)
    RZ(theta=-3.1655466193, wire=4)
  CNOT(control=2, target=4)
    RX(theta=-0.0166245844, wire=0)
    RY(theta=-0.1048111736, wire=0)
    RZ(theta=-1.6194782133, wire=0)
    RX(theta=-0.3613478766, wire=1)
    RY(theta=+3.2241531683, wire=1)
    RZ(theta=+2.2526573332, wire=1)
    RX(theta=-0.2293042453, wire=2)
    RY(theta=-2.6351030613, wire=2)
    RZ(theta=+1.4831543614, wire=2)
    RX(theta=+2.8280807642, wire=3)
    RY(theta=+0.3791126440, wire=3)
    RZ(theta=+1.9346443902, wire=3)
    RX(theta=-2.1442037622, wire=4)
    RY(theta=+0.8143111808, wire=4)
    RZ(theta=-0.0004603627, wire=4)
    RX(theta=+1.6837009213, wire=0)
    RY(theta=+2.4338981302, wire=0)
    RZ(theta=+1.8861260775, wire=0)
    RX(theta=-0.0111570389, wire=1)
    RY(theta=-3.0496351673, wire=1)
    RZ(theta=-0.7101684115, wire=1)
    RX(theta=-1.7363871421, wire=2)
    RY(theta=+2.1348825240, wire=2)
    RZ(theta=-3.0391414736, wire=2)
    RX(theta=-1.3268723351, wire=3)
    RY(theta=-2.6945484619, wire=3)
    RZ(theta=+0.9434351084, wire=3)
    RX(theta=+2.6412261708, wire=4)
    RY(theta=+2.4638982556, wire=4)
    RZ(theta=-2.1440574163, wire=4)
  CNOT(control=0, target=3)
    RX(theta=-0.0747853332, wire=0)
    RY(theta=-1.1560150274, wire=0)
    RZ(theta=-1.2525240229, wire=0)
    RX(theta=-0.5305170413, wire=1)
    RY(theta=-2.4660146624, wire=1)
    RZ(theta=-0.5946507258, wire=1)
    RX(theta=-1.0270448266, wire=2)
    RY(theta=+0.8914694855, wire=2)
    RZ(theta=+1.8465117211, wire=2)
    RX(theta=+0.9523948826, wire=3)
    RY(theta=+0.5511123103, wire=3)
    RZ(theta=-2.5498922914, wire=3)
    RX(theta=+0.9572139686, wire=4)
    RY(theta=-1.2729537054, wire=4)
    RZ(theta=+2.7374174845, wire=4)
  CNOT(control=1, target=3)
    RX(theta=-1.1005566390, wire=0)
    RY(theta=-2.0120710273, wire=0)
    RZ(theta=-0.9258176273, wire=0)
    RX(theta=-0.6694904324, wire=1)
    RY(theta=+2.4331594035, wire=1)
    RZ(theta=+1.0694563774, wire=1)
    RX(theta=-0.7570597612, wire=2)
    RY(theta=-2.1511025972, wire=2)
    RZ(theta=-2.7193153810, wire=2)
    RX(theta=-0.5405982463, wire=3)
    RY(theta=-0.6667029964, wire=3)
    RZ(theta=+2.2071360116, wire=3)
    RX(theta=+1.2578510489, wire=4)
    RY(theta=+2.6737534432, wire=4)
    RZ(theta=+1.5021394999, wire=4)
    RX(theta=+2.4353216988, wire=0)
    RY(theta=+1.4824929759, wire=0)
    RZ(theta=+2.8533946087, wire=0)
    RX(theta=-1.4317669721, wire=1)
    RY(theta=-1.8163200668, wire=1)
    RZ(theta=-2.8058735034, wire=1)
    RX(theta=+1.3244646679, wire=2)
    RY(theta=+1.0740898297, wire=2)
    RZ(theta=+0.4103483432, wire=2)
    RX(theta=+0.6406319286, wire=3)
    RY(theta=-1.4634778046, wire=3)
    RZ(theta=+1.2087323100, wire=3)
    RX(theta=-1.9452579850, wire=4)
    RY(theta=-0.4591595280, wire=4)
    RZ(theta=-0.0864790687, wire=4)
  CNOT(control=0, target=4)
    RX(theta=-2.1054542889, wire=0)
    RY(theta=+3.3735479828, wire=0)
    RZ(theta=-1.3514385786, wire=0)
    RX(theta=+2.2276149451, wire=1)
    RY(theta=-1.4645305732, wire=1)
    RZ(theta=-2.7354504464, wire=1)
    RX(theta=-0.7810461293, wire=2)
    RY(theta=-0.9278601130, wire=2)
    RZ(theta=+1.6543939208, wire=2)
    RX(theta=+0.1076926344, wire=3)
    RY(theta=+1.0224758273, wire=3)
    RZ(theta=+2.5947499836, wire=3)
    RX(theta=-1.6641500162, wire=4)
    RY(theta=+0.3884417421, wire=4)
    RZ(theta=-0.7012652496, wire=4)
  CNOT(control=2, target=4)
    RX(theta=+0.1459728471, wire=0)
    RY(theta=-1.0457753857, wire=0)
    RZ(theta=+1.0477227519, wire=0)
    RX(theta=+2.7901964530, wire=1)
    RY(theta=+2.8864337743, wire=1)
    RZ(theta=+3.0875589796, wire=1)
    RX(theta=-1.9193634611, wire=2)
    RY(theta=-1.0108562327, wire=2)
    RZ(theta=-0.8040589104, wire=2)
    RX(theta=+2.9792395546, wire=3)
    RY(theta=+2.4801727336, wire=3)
    RZ(theta=+0.1843776636, wire=3)
    RX(theta=+1.6066179238, wire=4)
    RY(theta=+0.2310811194, wire=4)
    RZ(theta=+1.1548435615, wire=4)
  CNOT(control=0, target=1)
    RX(theta=-3.1415926513, wire=0)
    RY(theta=+3.1415926505, wire=0)
    RZ(theta=+1.8472348907, wire=0)
    RX(theta=+1.8777871870, wire=1)
    RY(theta=+2.4320194496, wire=1)
    RZ(theta=+1.7261053842, wire=1)
    RX(theta=+0.6608051945, wire=2)
    RY(theta=+2.1214117264, wire=2)
    RZ(theta=-2.1967013966, wire=2)
    RX(theta=+2.7528226019, wire=3)
    RY(theta=+1.8455620876, wire=3)
    RZ(theta=-0.4891520793, wire=3)
    RX(theta=-2.1902454183, wire=4)
    RY(theta=+3.1156739568, wire=4)
    RZ(theta=+0.2239082604, wire=4)
  CNOT(control=0, target=3)
    RX(theta=+2.2069020315, wire=0)
    RY(theta=-1.7673471020, wire=0)
    RZ(theta=-0.4596529765, wire=0)
    RX(theta=+1.6432251151, wire=1)
    RY(theta=+0.0732466366, wire=1)
    RZ(theta=+1.6395011571, wire=1)
    RX(theta=+2.1748704387, wire=2)
    RY(theta=+3.4056801522, wire=2)
    RZ(theta=+0.6607016718, wire=2)
    RX(theta=+2.5832617027, wire=3)
    RY(theta=-2.6979426672, wire=3)
    RZ(theta=-2.4701073189, wire=3)
    RX(theta=-0.5731870551, wire=4)
    RY(theta=+2.8084504096, wire=4)
    RZ(theta=+2.4662902527, wire=4)
  CNOT(control=1, target=3)
    RX(theta=+1.6887013146, wire=0)
    RY(theta=-2.9755431231, wire=0)
    RZ(theta=+1.4446670103, wire=0)
    RX(theta=+2.0858547181, wire=1)
    RY(theta=-1.3569456887, wire=1)
    RZ(theta=+1.7584108219, wire=1)
    RX(theta=-2.2825667519, wire=2)
    RY(theta=-2.8751829497, wire=2)
    RZ(theta=+1.8996983500, wire=2)
    RX(theta=-3.1040828981, wire=3)
    RY(theta=+2.9003514453, wire=3)
    RZ(theta=-1.2078371088, wire=3)
    RX(theta=-2.3108168433, wire=4)
    RY(theta=+2.2836053549, wire=4)
    RZ(theta=-0.6390239106, wire=4)
  CNOT(control=0, target=2)
    RX(theta=+3.1415926552, wire=0)
    RY(theta=-3.1415926560, wire=0)
    RZ(theta=+1.4097015795, wire=0)
    RX(theta=+1.3591847209, wire=1)
    RY(theta=+0.8725463013, wire=1)
    RZ(theta=+1.8280212670, wire=1)
    RX(theta=+1.3555533196, wire=2)
    RY(theta=-2.6687556893, wire=2)
    RZ(theta=-0.6702572775, wire=2)
    RX(theta=-2.0550025934, wire=3)
    RY(theta=+1.9839438610, wire=3)
    RZ(theta=+1.0277026177, wire=3)
    RX(theta=+1.9583738052, wire=4)
    RY(theta=+1.6886989209, wire=4)
    RZ(theta=+2.8513027447, wire=4)
  CNOT(control=0, target=4)
    RX(theta=+1.3446661372, wire=0)
    RY(theta=+1.3335022225, wire=0)
    RZ(theta=+3.0149819139, wire=0)
    RX(theta=-0.5637167838, wire=1)
    RY(theta=-1.6982229809, wire=1)
    RZ(theta=+1.6144443566, wire=1)
    RX(theta=+2.1589595293, wire=2)
    RY(theta=-1.6204108769, wire=2)
    RZ(theta=+0.9401310608, wire=2)
    RX(theta=+2.2274030692, wire=3)
    RY(theta=+0.1088846953, wire=3)
    RZ(theta=+2.6626660368, wire=3)
    RX(theta=-2.8064584408, wire=4)
    RY(theta=+0.4438976585, wire=4)
    RZ(theta=-2.4702555277, wire=4)
  CNOT(control=2, target=4)
    RX(theta=+0.1575004631, wire=0)
    RY(theta=+1.9917358291, wire=0)
    RZ(theta=+1.8007132798, wire=0)
    RX(theta=+1.5814811376, wire=1)
    RY(theta=+2.8500813060, wire=1)
    RZ(theta=-2.1325160186, wire=1)
    RX(theta=+1.6142095506, wire=2)
    RY(theta=+1.5544354161, wire=2)
    RZ(theta=+0.9653736058, wire=2)
    RX(theta=+0.4747790713, wire=3)
    RY(theta=+3.0218485988, wire=3)
    RZ(theta=-2.4970208227, wire=3)
    RX(theta=+0.9169840744, wire=4)
    RY(theta=-0.8608532540, wire=4)
    RZ(theta=-2.2188003702, wire=4)
```

### (2) Merged form: net SU(2) per wire per block, alternating with CNOTs

**Block 0** (before CNOT (0, 1)):
```
  wire 0 (a (anc)): Rot(phi=-3.141593, theta=+1.570796, omega=+3.095612)
      U = [[+0.706920+0.016255j, +0.706920+0.016255j],
       [-0.706920+0.016255j, +0.706920-0.016255j]]
  wire 1 (     A1): Rot(phi=-0.192817, theta=+1.011321, omega=+2.124939)
      U = [[+0.497394-0.719701j, -0.193932+0.443869j],
       [+0.193932+0.443869j, +0.497394+0.719701j]]
  wire 2 (     A2): Rot(phi=-0.064497, theta=+1.442644, omega=-0.077567)
      U = [[+0.749040+0.053295j, -0.660364-0.004316j],
       [+0.660364-0.004316j, +0.749040-0.053295j]]
  wire 3 (     B1): Rot(phi=-0.194573, theta=+1.955564, omega=+0.407775)
      U = [[+0.555691-0.059462j, -0.791934+0.245993j],
       [+0.791934+0.245993j, +0.555691+0.059462j]]
  wire 4 (     B2): Rot(phi=+2.629551, theta=+1.667834, omega=-0.812771)
      U = [[+0.413278-0.529866j, +0.110936-0.732213j],
       [-0.110936-0.732213j, +0.413278+0.529866j]]
```
then **CNOT(control=0, target=1)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 1].

**Block 1** (before CNOT (1, 3)):
```
  wire 0 (a (anc)): Rot(phi=+3.381089, theta=+1.331781, omega=+2.412296)
      U = [[-0.762903-0.190662j, -0.546690-0.287675j],
       [+0.546690-0.287675j, -0.762903+0.190662j]]
  wire 1 (     A1): Rot(phi=+1.570796, theta=+2.844181, omega=-2.391226)
      U = [[+0.135867+0.059087j, +0.394405-0.906914j],
       [-0.394405-0.906914j, +0.135867-0.059087j]]
  wire 2 (     A2): Rot(phi=+0.976844, theta=+1.039949, omega=+4.368668)
      U = [[-0.774188-0.392129j, +0.062003+0.492974j],
       [-0.062003+0.492974j, -0.774188+0.392129j]]
  wire 3 (     B1): Rot(phi=-1.187353, theta=+1.878087, omega=+2.105774)
      U = [[+0.529380-0.261760j, +0.061085+0.804678j],
       [-0.061085+0.804678j, +0.529380+0.261760j]]
  wire 4 (     B2): Rot(phi=+1.998256, theta=+1.953478, omega=-1.574390)
      U = [[+0.547204-0.117739j, +0.177222-0.809504j],
       [-0.177222-0.809504j, +0.547204+0.117739j]]
```
then **CNOT(control=1, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [1, 3].

**Block 2** (before CNOT (0, 4)):
```
  wire 0 (a (anc)): Rot(phi=+3.870890, theta=+1.809812, omega=+0.917488)
      U = [[-0.453100-0.419914j, -0.073884-0.782888j],
       [+0.073884-0.782888j, -0.453100+0.419914j]]
  wire 1 (     A1): Rot(phi=+3.399424, theta=+1.569011, omega=-1.310465)
      U = [[+0.355534-0.611955j, +0.498928-0.500177j],
       [-0.498928-0.500177j, +0.355534+0.611955j]]
  wire 2 (     A2): Rot(phi=+1.575275, theta=+2.033424, omega=-2.010830)
      U = [[+0.513737+0.113683j, +0.187450-0.829465j],
       [-0.187450-0.829465j, +0.513737-0.113683j]]
  wire 3 (     B1): Rot(phi=+0.088635, theta=+2.043280, omega=+5.788665)
      U = [[-0.511256-0.105204j, +0.816963+0.245196j],
       [-0.816963+0.245196j, -0.511256+0.105204j]]
  wire 4 (     B2): Rot(phi=-0.310866, theta=+1.334705, omega=+0.869825)
      U = [[+0.754987-0.216674j, -0.514157+0.344514j],
       [+0.514157+0.344514j, +0.754987+0.216674j]]
```
then **CNOT(control=0, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 4].

**Block 3** (before CNOT (2, 4)):
```
  wire 0 (a (anc)): Rot(phi=-3.360873, theta=+1.563546, omega=+2.482307)
      U = [[+0.642288+0.301814j, +0.687557+0.153753j],
       [-0.687557+0.153753j, +0.642288-0.301814j]]
  wire 1 (     A1): Rot(phi=-3.074732, theta=+1.478454, omega=+2.657236)
      U = [[+0.722947+0.153145j, +0.648291+0.183341j],
       [-0.648291+0.183341j, +0.722947-0.153145j]]
  wire 2 (     A2): Rot(phi=+0.272631, theta=+1.689957, omega=+2.524559)
      U = [[+0.113734-0.653930j, -0.321851+0.675168j],
       [+0.321851+0.675168j, +0.113734+0.653930j]]
  wire 3 (     B1): Rot(phi=-2.251779, theta=+1.553682, omega=+0.115308)
      U = [[+0.343495+0.624954j, -0.264742+0.649119j],
       [+0.264742+0.649119j, +0.343495-0.624954j]]
  wire 4 (     B2): Rot(phi=+0.000000, theta=+1.570796, omega=-2.536521)
      U = [[+0.210677+0.674993j, -0.210677-0.674993j],
       [+0.210677-0.674993j, +0.210677-0.674993j]]
```
then **CNOT(control=2, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [2, 4].

**Block 4** (before CNOT (0, 3)):
```
  wire 0 (a (anc)): Rot(phi=-2.482307, theta=+1.578046, omega=+0.292979)
      U = [[+0.322922+0.626176j, -0.129252+0.697796j],
       [+0.129252+0.697796j, +0.322922-0.626176j]]
  wire 1 (     A1): Rot(phi=-1.493318, theta=+0.440927, omega=-1.467170)
      U = [[+0.088240+0.971798j, -0.218663+0.002859j],
       [+0.218663+0.002859j, +0.088240-0.971798j]]
  wire 2 (     A2): Rot(phi=+1.209761, theta=+1.657328, omega=-3.915633)
      U = [[+0.146081+0.659885j, +0.616943-0.403228j],
       [-0.616943-0.403228j, +0.146081-0.659885j]]
  wire 3 (     B1): Rot(phi=+2.358674, theta=+1.390469, omega=-1.015773)
      U = [[+0.601207-0.477730j, +0.074411-0.636229j],
       [-0.074411-0.636229j, +0.601207+0.477730j]]
  wire 4 (     B2): Rot(phi=+0.203900, theta=+1.898010, omega=-2.774986)
      U = [[+0.163914+0.558954j, -0.066054-0.810148j],
       [+0.066054-0.810148j, +0.163914-0.558954j]]
```
then **CNOT(control=0, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 3].

**Block 5** (before CNOT (1, 3)):
```
  wire 0 (a (anc)): Rot(phi=-3.108705, theta=+1.157245, omega=+1.807384)
      U = [[+0.666160+0.507112j, +0.424009+0.345376j],
       [-0.424009+0.345376j, +0.666160-0.507112j]]
  wire 1 (     A1): Rot(phi=-3.704774, theta=+2.309167, omega=+1.793495)
      U = [[+0.233316+0.330185j, +0.845089+0.349809j],
       [-0.845089+0.349809j, +0.233316-0.330185j]]
  wire 2 (     A2): Rot(phi=-0.604714, theta=+1.239747, omega=+2.977677)
      U = [[+0.305170-0.754579j, +0.127003+0.566880j],
       [-0.127003+0.566880j, +0.305170+0.754579j]]
  wire 3 (     B1): Rot(phi=+0.924527, theta=+1.054228, omega=-3.764054)
      U = [[+0.130036+0.854424j, +0.351444-0.359913j],
       [-0.351444-0.359913j, +0.130036-0.854424j]]
  wire 4 (     B2): Rot(phi=+2.895693, theta=+1.401008, omega=+0.574070)
      U = [[-0.124884-0.754248j, -0.256936-0.591183j],
       [+0.256936-0.591183j, -0.124884+0.754248j]]
```
then **CNOT(control=1, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [1, 3].

**Block 6** (before CNOT (0, 4)):
```
  wire 0 (a (anc)): Rot(phi=-1.807384, theta=+1.984347, omega=+2.234522)
      U = [[+0.534447-0.115909j, +0.364278+0.753812j],
       [-0.364278+0.753812j, +0.534447+0.115909j]]
  wire 1 (     A1): Rot(phi=-1.053618, theta=+1.661416, omega=-3.316109)
      U = [[-0.388560+0.551156j, -0.314217-0.668218j],
       [+0.314217-0.668218j, -0.388560-0.551156j]]
  wire 2 (     A2): Rot(phi=-0.626528, theta=+2.842543, omega=+2.011019)
      U = [[+0.114678-0.095082j, -0.246581+0.957604j],
       [+0.246581+0.957604j, +0.114678+0.095082j]]
  wire 3 (     B1): Rot(phi=+5.093734, theta=+1.658195, omega=-0.475251)
      U = [[-0.454733-0.499574j, +0.690807-0.257737j],
       [-0.690807-0.257737j, -0.454733+0.499574j]]
  wire 4 (     B2): Rot(phi=-2.892834, theta=+1.099666, omega=-0.409118)
      U = [[-0.068289+0.849873j, -0.168802+0.494529j],
       [+0.168802+0.494529j, -0.068289-0.849873j]]
```
then **CNOT(control=0, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 4].

**Block 7** (before CNOT (2, 4)):
```
  wire 0 (a (anc)): Rot(phi=+1.838714, theta=+1.051925, omega=+0.084056)
      U = [[+0.495022-0.709157j, -0.320911-0.386090j],
       [+0.320911-0.386090j, +0.495022+0.709157j]]
  wire 1 (     A1): Rot(phi=-3.225867, theta=+1.635606, omega=+2.631024)
      U = [[+0.653803+0.200399j, +0.713133+0.154347j],
       [-0.713133+0.154347j, +0.653803-0.200399j]]
  wire 2 (     A2): Rot(phi=+3.626916, theta=+1.130967, omega=-2.378797)
      U = [[+0.685185-0.493371j, +0.530676-0.074100j],
       [-0.530676-0.074100j, +0.685185+0.493371j]]
  wire 3 (     B1): Rot(phi=+0.065557, theta=+1.026011, omega=+2.468739)
      U = [[+0.260514-0.831414j, -0.177116+0.457725j],
       [+0.177116+0.457725j, +0.260514+0.831414j]]
  wire 4 (     B2): Rot(phi=-1.180823, theta=+1.657177, omega=+0.904977)
      U = [[+0.669498+0.092929j, -0.371229+0.636652j],
       [+0.371229+0.636652j, +0.669498-0.092929j]]
```
then **CNOT(control=2, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [2, 4].

**Block 8** (before CNOT (0, 1)):
```
  wire 0 (a (anc)): Rot(phi=+3.057537, theta=+1.051925, omega=-1.925575)
      U = [[+0.729980-0.463766j, +0.399659-0.303846j],
       [-0.399659-0.303846j, +0.729980+0.463766j]]
  wire 1 (     A1): Rot(phi=-0.922317, theta=+0.431129, omega=+0.913820)
      U = [[+0.976847+0.004150j, -0.129913+0.169928j],
       [+0.129913+0.169928j, +0.976847-0.004150j]]
  wire 2 (     A2): Rot(phi=-2.609176, theta=+1.753216, omega=+0.468041)
      U = [[+0.306847+0.561373j, -0.024735+0.768175j],
       [+0.024735+0.768175j, +0.306847-0.561373j]]
  wire 3 (     B1): Rot(phi=-0.204752, theta=+0.678136, omega=-2.696616)
      U = [[+0.113002+0.936271j, -0.106162-0.315211j],
       [+0.106162-0.315211j, +0.113002-0.936271j]]
  wire 4 (     B2): Rot(phi=+1.339572, theta=+1.605665, omega=-0.424160)
      U = [[+0.623168-0.306969j, -0.457284-0.555267j],
       [+0.457284-0.555267j, +0.623168+0.306969j]]
```
then **CNOT(control=0, target=1)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 1].

**Block 9** (before CNOT (0, 3)):
```
  wire 0 (a (anc)): Rot(phi=+0.637331, theta=+0.000000, omega=+4.351497)
      U = [[-0.797788-0.602938j, +0.000000+0.000000j],
       [-0.000000+0.000000j, -0.797788+0.602938j]]
  wire 1 (     A1): Rot(phi=-0.837480, theta=+1.339485, omega=-0.048364)
      U = [[+0.708330+0.336000j, -0.573087+0.238630j],
       [+0.573087+0.238630j, +0.708330-0.336000j]]
  wire 2 (     A2): Rot(phi=-0.360361, theta=+1.996624, omega=-2.936227)
      U = [[-0.041940+0.540096j, -0.234605-0.807154j],
       [+0.234605-0.807154j, -0.041940-0.540096j]]
  wire 3 (     B1): Rot(phi=-0.106449, theta=+1.317006, omega=-3.228381)
      U = [[-0.076298+0.787220j, -0.006015-0.611904j],
       [+0.006015-0.611904j, -0.076298-0.787220j]]
  wire 4 (     B2): Rot(phi=+1.538967, theta=+0.951587, omega=+1.813182)
      U = [[-0.093412-0.884008j, -0.453745+0.062605j],
       [+0.453745+0.062605j, -0.093412+0.884008j]]
```
then **CNOT(control=0, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 3].

**Block 10** (before CNOT (1, 3)):
```
  wire 0 (a (anc)): Rot(phi=+3.300420, theta=+1.454520, omega=-1.403607)
      U = [[+0.435484-0.606927j, +0.468133-0.472064j],
       [-0.468133-0.472064j, +0.435484+0.606927j]]
  wire 1 (     A1): Rot(phi=+1.497358, theta=+1.643031, omega=+0.063395)
      U = [[+0.484032-0.479195j, -0.551911-0.481124j],
       [+0.551911-0.481124j, +0.484032+0.479195j]]
  wire 2 (     A2): Rot(phi=-1.888233, theta=+0.990456, omega=-0.731861)
      U = [[+0.226832+0.850118j, -0.397986+0.259717j],
       [+0.397986+0.259717j, +0.226832-0.850118j]]
  wire 3 (     B1): Rot(phi=-2.301996, theta=+0.698150, omega=+2.844360)
      U = [[+0.905348-0.251715j, +0.288246+0.184113j],
       [-0.288246+0.184113j, +0.905348+0.251715j]]
  wire 4 (     B2): Rot(phi=+1.002845, theta=+2.488130, omega=+3.568152)
      U = [[-0.210348-0.242409j, -0.269139+0.908051j],
       [+0.269139+0.908051j, -0.210348+0.242409j]]
```
then **CNOT(control=1, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [1, 3].

**Block 11** (before CNOT (0, 2)):
```
  wire 0 (a (anc)): Rot(phi=+4.545200, theta=+1.454520, omega=-0.106553)
      U = [[-0.451197-0.595339j, +0.455636-0.484137j],
       [-0.455636-0.484137j, -0.451197+0.595339j]]
  wire 1 (     A1): Rot(phi=+2.954800, theta=+1.675526, omega=+0.692854)
      U = [[-0.167509-0.647821j, -0.316418-0.672420j],
       [+0.316418-0.672420j, -0.167509+0.647821j]]
  wire 2 (     A2): Rot(phi=+1.916705, theta=+0.889072, omega=-3.036002)
      U = [[+0.765079+0.479292j, +0.338342-0.265440j],
       [-0.338342-0.265440j, +0.765079-0.479292j]]
  wire 3 (     B1): Rot(phi=+0.151258, theta=+0.244084, omega=+1.777950)
      U = [[+0.565504-0.815711j, -0.083643+0.088454j],
       [+0.083643+0.088454j, +0.565504+0.815711j]]
  wire 4 (     B2): Rot(phi=+0.568134, theta=+1.114119, omega=+1.536301)
      U = [[+0.420711-0.737215j, -0.467947+0.246052j],
       [+0.467947+0.246052j, +0.420711+0.737215j]]
```
then **CNOT(control=0, target=2)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 2].

**Block 12** (before CNOT (0, 4)):
```
  wire 0 (a (anc)): Rot(phi=+0.602299, theta=+0.000000, omega=+3.948995)
      U = [[-0.647920-0.761708j, +0.000000+0.000000j],
       [-0.000000+0.000000j, -0.647920+0.761708j]]
  wire 1 (     A1): Rot(phi=+0.687164, theta=+1.435355, omega=+0.420315)
      U = [[+0.640760-0.396157j, -0.651792-0.087485j],
       [+0.651792-0.087485j, +0.640760+0.396157j]]
  wire 2 (     A2): Rot(phi=-2.053163, theta=+1.762112, omega=+3.942892)
      U = [[+0.372800-0.515699j, +0.763475+0.110368j],
       [-0.763475+0.110368j, +0.372800+0.515699j]]
  wire 3 (     B1): Rot(phi=+0.370106, theta=+1.382793, omega=+3.047412)
      U = [[-0.105943-0.763037j, -0.146692+0.620510j],
       [+0.146692+0.620510j, -0.105943+0.763037j]]
  wire 4 (     B2): Rot(phi=-0.109230, theta=+1.526324, omega=+0.895361)
      U = [[+0.667544-0.276792j, -0.605832+0.332775j],
       [+0.605832+0.332775j, +0.667544+0.276792j]]
```
then **CNOT(control=0, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 4].

**Block 13** (before CNOT (2, 4)):
```
  wire 0 (a (anc)): Rot(phi=+0.231469, theta=+1.518067, omega=+1.664184)
      U = [[+0.423294-0.589216j, -0.519058+0.451914j],
       [+0.519058+0.451914j, +0.423294+0.589216j]]
  wire 1 (     A1): Rot(phi=+3.073240, theta=+1.678423, omega=-2.094548)
      U = [[+0.589648-0.314016j, +0.631365-0.393812j],
       [-0.631365-0.393812j, +0.589648+0.314016j]]
  wire 2 (     A2): Rot(phi=+3.182880, theta=+1.543276, omega=-0.043070)
      U = [[+0.000639-0.716769j, +0.029403-0.696690j],
       [-0.029403-0.696690j, +0.000639+0.716769j]]
  wire 3 (     B1): Rot(phi=+1.433648, theta=+2.222847, omega=+1.008315)
      U = [[+0.151959-0.416533j, -0.876138-0.189186j],
       [+0.876138-0.189186j, +0.151959+0.416533j]]
  wire 4 (     B2): Rot(phi=-0.605073, theta=+2.592201, omega=-0.010043)
      U = [[+0.258526+0.082117j, -0.920223+0.282154j],
       [+0.920223+0.282154j, +0.258526-0.082117j]]
```
then **CNOT(control=2, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [2, 4].

**Block 14** (before measurement):
```
  wire 0 (a (anc)): Rot(phi=-0.070107, theta=+1.986201, omega=+1.628433)
      U = [[+0.388548-0.383732j, -0.553343+0.628961j],
       [+0.553343+0.628961j, +0.388548+0.383732j]]
  wire 1 (     A1): Rot(phi=-1.279269, theta=+1.560562, omega=-3.706383)
      U = [[-0.566319+0.429417j, -0.245999-0.659066j],
       [+0.245999-0.659066j, -0.566319-0.429417j]]
  wire 2 (     A2): Rot(phi=+0.016345, theta=+1.571506, omega=-0.648830)
      U = [[+0.671803+0.219830j, -0.668595-0.230945j],
       [+0.668595-0.230945j, +0.671803-0.219830j]]
  wire 3 (     B1): Rot(phi=-1.313433, theta=+2.653065, omega=-3.839459)
      U = [[-0.204238+0.129516j, -0.293954-0.924718j],
       [+0.293954-0.924718j, -0.204238-0.129516j]]
  wire 4 (     B2): Rot(phi=-3.740266, theta=+1.163171, omega=+1.967172)
      U = [[+0.528171+0.647495j, +0.526744+0.155968j],
       [-0.526744+0.155968j, +0.528171-0.647495j]]
```

==============================================================================
## 5c — 5-CNOT observable circuit
==============================================================================

Pruned from the 14-CNOT circuit under the OBSERVABLE cost: it need only reproduce the ancilla-parity correlators <Z_a> -> Tr(rho^2) and <Z_a (x) O> -> Tr(O rho^2) for O in {|Phi+><Phi+|, ZZ} over an eps grid (NOT the full unitary). This is the 5c rung of the relaxation ladder; 4 CNOTs is unreachable.

- CNOT sequence (in order): [(0, 1), (2, 4), (0, 4), (2, 4), (0, 3)]
- CNOT count: 5    rotation blocks: 6
- verification: |Tr(U_target^dag V)|/32 = 0.056068843033  (NOT ~1: it is not the full unitary, by design)
- observable read-out F = <Z_a (x) O> / <Z_a>, with the two register copies measured jointly; O acts on kept register A = wires (1,2).
- drawn(merged)==solution overlap = 1.000000000000

### (1) Exact primitive gate list
```
    RX(theta=-0.2452369203, wire=0)
    RY(theta=-2.1122366701, wire=0)
    RZ(theta=-1.7817128909, wire=0)
    RX(theta=-0.5402475723, wire=1)
    RY(theta=+1.1359438838, wire=1)
    RZ(theta=+1.9465485012, wire=1)
    RX(theta=+2.4345111147, wire=2)
    RY(theta=+1.8004517209, wire=2)
    RZ(theta=+2.8584166226, wire=2)
    RX(theta=+0.7052703324, wire=3)
    RY(theta=+2.2232951458, wire=3)
    RZ(theta=+0.7053206136, wire=3)
    RX(theta=+2.0267365980, wire=4)
    RY(theta=-0.6390442243, wire=4)
    RZ(theta=+0.9440466200, wire=4)
  CNOT(control=0, target=1)
    RX(theta=+1.6065574747, wire=0)
    RY(theta=-0.3777804278, wire=0)
    RZ(theta=+2.9814829684, wire=0)
    RX(theta=+0.8167076994, wire=1)
    RY(theta=+0.8581881053, wire=1)
    RZ(theta=-2.3910801210, wire=1)
    RX(theta=-0.6803909881, wire=2)
    RY(theta=+1.3346373934, wire=2)
    RZ(theta=+1.5968334963, wire=2)
    RX(theta=-2.2602406599, wire=3)
    RY(theta=-1.7839391087, wire=3)
    RZ(theta=+1.1455998535, wire=3)
    RX(theta=+1.3058131542, wire=4)
    RY(theta=-0.4631145515, wire=4)
    RZ(theta=+1.5715059747, wire=4)
    RX(theta=-2.5065259933, wire=0)
    RY(theta=+2.2034137133, wire=0)
    RZ(theta=-1.6279128402, wire=0)
    RX(theta=+1.3288477380, wire=1)
    RY(theta=-2.4461534437, wire=1)
    RZ(theta=-1.2284003888, wire=1)
    RX(theta=-0.3866928950, wire=2)
    RY(theta=+1.2750654401, wire=2)
    RZ(theta=+2.8243188025, wire=2)
    RX(theta=-1.6640056193, wire=3)
    RY(theta=+0.3031556098, wire=3)
    RZ(theta=-1.1715755631, wire=3)
    RX(theta=+0.7143480129, wire=4)
    RY(theta=-0.0258434855, wire=4)
    RZ(theta=-2.1860893018, wire=4)
    RX(theta=+2.1187395943, wire=0)
    RY(theta=-1.2486743311, wire=0)
    RZ(theta=+3.4328267094, wire=0)
    RX(theta=-2.9989068318, wire=1)
    RY(theta=-0.0497948666, wire=1)
    RZ(theta=-2.5265327369, wire=1)
    RX(theta=-0.7678584673, wire=2)
    RY(theta=-0.8970089399, wire=2)
    RZ(theta=+1.4912670853, wire=2)
    RX(theta=-0.3293541784, wire=3)
    RY(theta=+2.5321022887, wire=3)
    RZ(theta=+2.4953728804, wire=3)
    RX(theta=-0.2210470685, wire=4)
    RY(theta=+1.7461201198, wire=4)
    RZ(theta=-2.4826179865, wire=4)
    RX(theta=-0.7980975463, wire=0)
    RY(theta=+0.4580095976, wire=0)
    RZ(theta=-1.4287241091, wire=0)
    RX(theta=+2.5577503334, wire=1)
    RY(theta=+3.2401257940, wire=1)
    RZ(theta=-2.7669994583, wire=1)
    RX(theta=+2.0790555989, wire=2)
    RY(theta=+0.2115785616, wire=2)
    RZ(theta=+0.8312351986, wire=2)
    RX(theta=-3.1282125703, wire=3)
    RY(theta=+3.2230145731, wire=3)
    RZ(theta=-0.7946364969, wire=3)
    RX(theta=-2.6743318615, wire=4)
    RY(theta=+2.0112592239, wire=4)
    RZ(theta=+1.3481039560, wire=4)
    RX(theta=-1.6654925941, wire=0)
    RY(theta=-1.6266250649, wire=0)
    RZ(theta=+0.9367435522, wire=0)
    RX(theta=-0.1209658470, wire=1)
    RY(theta=-1.8023111540, wire=1)
    RZ(theta=+0.6412958161, wire=1)
    RX(theta=-1.2324983291, wire=2)
    RY(theta=+1.6543156948, wire=2)
    RZ(theta=+1.4267407966, wire=2)
    RX(theta=-1.1853138522, wire=3)
    RY(theta=-0.8956964292, wire=3)
    RZ(theta=-1.1581979994, wire=3)
    RX(theta=-0.4186109011, wire=4)
    RY(theta=+1.5725275323, wire=4)
    RZ(theta=-3.4577427346, wire=4)
  CNOT(control=2, target=4)
    RX(theta=-0.1296269549, wire=0)
    RY(theta=-0.0647691819, wire=0)
    RZ(theta=-1.5612615359, wire=0)
    RX(theta=-0.4187486030, wire=1)
    RY(theta=+3.0672231279, wire=1)
    RZ(theta=+1.6715072312, wire=1)
    RX(theta=-0.5223448549, wire=2)
    RY(theta=-2.6318658561, wire=2)
    RZ(theta=+1.3319295500, wire=2)
    RX(theta=+2.6730544359, wire=3)
    RY(theta=+0.4731391770, wire=3)
    RZ(theta=+1.6382485248, wire=3)
    RX(theta=-1.9095709525, wire=4)
    RY(theta=+1.5336007487, wire=4)
    RZ(theta=+0.1240100510, wire=4)
    RX(theta=+1.6754689717, wire=0)
    RY(theta=+2.2752419945, wire=0)
    RZ(theta=+1.8135281602, wire=0)
    RX(theta=+0.0107047132, wire=1)
    RY(theta=-2.7886941434, wire=1)
    RZ(theta=-0.2065307534, wire=1)
    RX(theta=-1.6829957326, wire=2)
    RY(theta=+1.9743276577, wire=2)
    RZ(theta=-2.9796772806, wire=2)
    RX(theta=-1.3072810475, wire=3)
    RY(theta=-3.0481390034, wire=3)
    RZ(theta=+0.7341159587, wire=3)
    RX(theta=+3.0877919660, wire=4)
    RY(theta=+1.8165294141, wire=4)
    RZ(theta=-2.3623328851, wire=4)
    RX(theta=-0.1015306618, wire=0)
    RY(theta=-1.1563557496, wire=0)
    RZ(theta=-1.2988258484, wire=0)
    RX(theta=-0.1831006129, wire=1)
    RY(theta=-2.1439307846, wire=1)
    RZ(theta=-0.4660147519, wire=1)
    RX(theta=-0.7720700625, wire=2)
    RY(theta=+1.0057291775, wire=2)
    RZ(theta=+1.5913903198, wire=2)
    RX(theta=+1.1341444462, wire=3)
    RY(theta=-0.0022907136, wire=3)
    RZ(theta=-3.4345906698, wire=3)
    RX(theta=+0.7682600527, wire=4)
    RY(theta=-0.7469571969, wire=4)
    RZ(theta=+2.9876764164, wire=4)
    RX(theta=-1.1015141320, wire=0)
    RY(theta=-2.0904466379, wire=0)
    RZ(theta=-0.9804719107, wire=0)
    RX(theta=-0.2717452930, wire=1)
    RY(theta=+2.3809684522, wire=1)
    RZ(theta=+0.7468204228, wire=1)
    RX(theta=-0.8786544342, wire=2)
    RY(theta=-2.2401783026, wire=2)
    RZ(theta=-2.6618607056, wire=2)
    RX(theta=-0.6158483660, wire=3)
    RY(theta=-0.3126130021, wire=3)
    RZ(theta=+2.3299720803, wire=3)
    RX(theta=+1.7834503828, wire=4)
    RY(theta=+2.1966496141, wire=4)
    RZ(theta=+1.3531395164, wire=4)
    RX(theta=+2.3322357747, wire=0)
    RY(theta=+1.5171591112, wire=0)
    RZ(theta=+2.9652858017, wire=0)
    RX(theta=-1.4601317902, wire=1)
    RY(theta=-2.1346529508, wire=1)
    RZ(theta=-2.9000581373, wire=1)
    RX(theta=+1.0408330136, wire=2)
    RY(theta=+0.9955306874, wire=2)
    RZ(theta=+0.6212663386, wire=2)
    RX(theta=+0.4932613947, wire=3)
    RY(theta=-1.8741928658, wire=3)
    RZ(theta=+1.0868827660, wire=3)
    RX(theta=-1.4203151144, wire=4)
    RY(theta=-0.6985398307, wire=4)
    RZ(theta=+0.6675458198, wire=4)
  CNOT(control=0, target=4)
    RX(theta=-2.0251262398, wire=0)
    RY(theta=+3.5302422356, wire=0)
    RZ(theta=-1.2197841395, wire=0)
    RX(theta=+2.5256961443, wire=1)
    RY(theta=-1.6716245250, wire=1)
    RZ(theta=-2.3827604869, wire=1)
    RX(theta=-0.9147118794, wire=2)
    RY(theta=-0.8959914135, wire=2)
    RZ(theta=+1.6765991312, wire=2)
    RX(theta=+0.5473331820, wire=3)
    RY(theta=+0.9957802011, wire=3)
    RZ(theta=+2.1396755526, wire=3)
    RX(theta=-1.2433493079, wire=4)
    RY(theta=-0.7117938707, wire=4)
    RZ(theta=-1.5763299737, wire=4)
  CNOT(control=2, target=4)
    RX(theta=+0.3213700723, wire=0)
    RY(theta=-1.1102943572, wire=0)
    RZ(theta=+1.3452473641, wire=0)
    RX(theta=+2.8413291636, wire=1)
    RY(theta=+2.5042788737, wire=1)
    RZ(theta=+3.4052724871, wire=1)
    RX(theta=-1.6219547685, wire=2)
    RY(theta=-0.8769650030, wire=2)
    RZ(theta=-0.2267744606, wire=2)
    RX(theta=+2.8857910484, wire=3)
    RY(theta=+2.3645180182, wire=3)
    RZ(theta=-0.1086681961, wire=3)
    RX(theta=+1.2376606417, wire=4)
    RY(theta=+0.2399836136, wire=4)
    RZ(theta=+1.2647091610, wire=4)
    RX(theta=-3.1260703628, wire=0)
    RY(theta=+3.1279282631, wire=0)
    RZ(theta=+1.7957508987, wire=0)
    RX(theta=+1.9499459768, wire=1)
    RY(theta=+2.4090388994, wire=1)
    RZ(theta=+1.8890236120, wire=1)
    RX(theta=+0.4870415179, wire=2)
    RY(theta=+1.9645980035, wire=2)
    RZ(theta=-2.3255580478, wire=2)
    RX(theta=+3.1529955298, wire=3)
    RY(theta=+2.0847974396, wire=3)
    RZ(theta=-0.9611265301, wire=3)
    RX(theta=-2.3874907450, wire=4)
    RY(theta=+3.4326609239, wire=4)
    RZ(theta=+0.0079633753, wire=4)
  CNOT(control=0, target=3)
    RX(theta=+2.1603555673, wire=0)
    RY(theta=-1.7607885212, wire=0)
    RZ(theta=-0.4615562473, wire=0)
    RX(theta=+1.5446214668, wire=1)
    RY(theta=-0.1376415177, wire=1)
    RZ(theta=+1.3699277555, wire=1)
    RX(theta=+1.6210686664, wire=2)
    RY(theta=+3.6715177172, wire=2)
    RZ(theta=+0.6652857781, wire=2)
    RX(theta=+2.8087206615, wire=3)
    RY(theta=-3.0794947806, wire=3)
    RZ(theta=-3.1023815251, wire=3)
    RX(theta=-0.3689567935, wire=4)
    RY(theta=+3.0366820259, wire=4)
    RZ(theta=+2.7402829772, wire=4)
    RX(theta=+1.6155319705, wire=0)
    RY(theta=-2.9874652690, wire=0)
    RZ(theta=+1.4250577284, wire=0)
    RX(theta=+1.9159056497, wire=1)
    RY(theta=-0.9308372077, wire=1)
    RZ(theta=+1.3655524280, wire=1)
    RX(theta=-1.9429284172, wire=2)
    RY(theta=-3.2285009263, wire=2)
    RZ(theta=+2.3531814743, wire=2)
    RX(theta=-3.1040828981, wire=3)
    RY(theta=+2.9003514453, wire=3)
    RZ(theta=-1.2078371088, wire=3)
    RX(theta=-2.1758088891, wire=4)
    RY(theta=+2.7141126996, wire=4)
    RZ(theta=-0.8748092795, wire=4)
    RX(theta=+3.1420115375, wire=0)
    RY(theta=-3.3155680712, wire=0)
    RZ(theta=+1.3624362601, wire=0)
    RX(theta=+0.9189048998, wire=1)
    RY(theta=+1.3223184445, wire=1)
    RZ(theta=+2.0127318739, wire=1)
    RX(theta=+1.8158986457, wire=2)
    RY(theta=-3.0995749707, wire=2)
    RZ(theta=-0.4723659584, wire=2)
    RX(theta=-2.0550025934, wire=3)
    RY(theta=+1.9839438610, wire=3)
    RZ(theta=+1.0277026177, wire=3)
    RX(theta=+2.2941180370, wire=4)
    RY(theta=+1.7873167236, wire=4)
    RZ(theta=+2.4346005287, wire=4)
    RX(theta=+1.0453855037, wire=0)
    RY(theta=+1.0640349406, wire=0)
    RZ(theta=+3.5341307577, wire=0)
    RX(theta=-0.8541603840, wire=1)
    RY(theta=-2.2546170436, wire=1)
    RZ(theta=+1.1054167561, wire=1)
    RX(theta=+1.5622262971, wire=2)
    RY(theta=-1.7615810267, wire=2)
    RZ(theta=+0.3668931423, wire=2)
    RX(theta=+2.2274030692, wire=3)
    RY(theta=+0.1088846953, wire=3)
    RZ(theta=+2.6626660368, wire=3)
    RX(theta=-3.0504846445, wire=4)
    RY(theta=+0.1738301895, wire=4)
    RZ(theta=-2.4273407793, wire=4)
    RX(theta=-0.2697003057, wire=0)
    RY(theta=+2.2237797354, wire=0)
    RZ(theta=+1.8007132798, wire=0)
    RX(theta=+2.1705894715, wire=1)
    RY(theta=+3.3952119965, wire=1)
    RZ(theta=-2.2469443660, wire=1)
    RX(theta=+1.7650829158, wire=2)
    RY(theta=+2.1342313420, wire=2)
    RZ(theta=+0.8509452583, wire=2)
    RX(theta=+0.4747790713, wire=3)
    RY(theta=+3.0218485988, wire=3)
    RZ(theta=-2.4970208227, wire=3)
    RX(theta=+0.9169840744, wire=4)
    RY(theta=-0.8608532540, wire=4)
    RZ(theta=-2.2188003702, wire=4)
```

### (2) Merged form: net SU(2) per wire per block, alternating with CNOTs

**Block 0** (before CNOT (0, 1)):
```
  wire 0 (a (anc)): Rot(phi=-3.286577, theta=+2.094338, omega=+1.075735)
      U = [[+0.224390+0.446849j, +0.496365+0.709646j],
       [-0.496365+0.709646j, +0.224390-0.446849j]]
  wire 1 (     A1): Rot(phi=-0.234522, theta=+1.201157, omega=+2.530836)
      U = [[+0.338393-0.752416j, -0.105683+0.555150j],
       [+0.105683+0.555150j, +0.338393+0.752416j]]
  wire 2 (     A2): Rot(phi=-0.150716, theta=+1.396853, omega=+0.437071)
      U = [[+0.758018-0.109279j, -0.615443+0.186269j],
       [+0.615443+0.186269j, +0.758018+0.109279j]]
  wire 3 (     B1): Rot(phi=-0.459925, theta=+2.051411, omega=-0.114555)
      U = [[+0.497253+0.146893j, -0.842363+0.146927j],
       [+0.842363+0.146927j, +0.497253-0.146893j]]
  wire 4 (     B2): Rot(phi=+2.262142, theta=+1.932020, omega=-0.342197)
      U = [[+0.326108-0.465772j, -0.218331-0.793121j],
       [+0.218331-0.793121j, +0.326108+0.465772j]]
```
then **CNOT(control=0, target=1)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 1].

**Block 1** (before CNOT (2, 4)):
```
  wire 0 (a (anc)): Rot(phi=-4.568901, theta=+1.571564, omega=-0.736983)
      U = [[-0.624112+0.331813j, +0.239341+0.665657j],
       [-0.239341+0.665657j, -0.624112-0.331813j]]
  wire 1 (     A1): Rot(phi=+1.334259, theta=+1.498071, omega=+3.273327)
      U = [[-0.490012-0.544260j, -0.385192+0.561512j],
       [+0.385192+0.561512j, -0.490012+0.544260j]]
  wire 2 (     A2): Rot(phi=+0.263546, theta=+1.644722, omega=+6.001296)
      U = [[-0.680464-0.006241j, +0.705674+0.197367j],
       [-0.705674+0.197367j, -0.680464+0.006241j]]
  wire 3 (     B1): Rot(phi=+5.078465, theta=+1.246326, omega=+0.152645)
      U = [[-0.702252-0.407732j, +0.454279-0.366370j],
       [-0.454279-0.366370j, -0.702252+0.407732j]]
  wire 4 (     B2): Rot(phi=+1.627770, theta=+1.391370, omega=-3.388649)
      U = [[+0.488827+0.591846j, +0.516594-0.379339j],
       [-0.516594-0.379339j, +0.488827-0.591846j]]
```
then **CNOT(control=2, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [2, 4].

**Block 2** (before CNOT (0, 4)):
```
  wire 0 (a (anc)): Rot(phi=+0.737252, theta=+1.569133, omega=-3.905372)
      U = [[-0.009386+0.707632j, +0.481857-0.516703j],
       [-0.481857-0.516703j, -0.009386-0.707632j]]
  wire 1 (     A1): Rot(phi=+2.035311, theta=+1.009143, omega=+4.240231)
      U = [[-0.875375-0.003345j, -0.218222+0.431377j],
       [+0.218222+0.431377j, -0.875375+0.003345j]]
  wire 2 (     A2): Rot(phi=-0.354518, theta=+1.962606, omega=+2.576398)
      U = [[+0.246737-0.498187j, -0.087398+0.826615j],
       [+0.087398+0.826615j, +0.246737+0.498187j]]
  wire 3 (     B1): Rot(phi=-0.917129, theta=+1.313598, omega=-4.364923)
      U = [[-0.694787+0.380075j, +0.093116-0.603443j],
       [-0.093116-0.603443j, -0.694787-0.380075j]]
  wire 4 (     B2): Rot(phi=-0.004938, theta=+1.246324, omega=-1.548920)
      U = [[+0.579040+0.569313j, -0.418167-0.407102j],
       [+0.418167-0.407102j, +0.579040-0.569313j]]
```
then **CNOT(control=0, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 4].

**Block 3** (before CNOT (2, 4)):
```
  wire 0 (a (anc)): Rot(phi=+1.998385, theta=+1.152580, omega=+0.168007)
      U = [[+0.392839-0.740772j, -0.332203-0.431944j],
       [+0.332203-0.431944j, +0.392839+0.740772j]]
  wire 1 (     A1): Rot(phi=+3.199972, theta=+1.488541, omega=-3.001060)
      U = [[+0.731947-0.073037j, +0.676864-0.027819j],
       [-0.676864-0.027819j, +0.731947+0.073037j]]
  wire 2 (     A2): Rot(phi=+3.706632, theta=+1.179804, omega=-2.494547)
      U = [[+0.683002-0.473351j, +0.555812-0.022803j],
       [-0.555812-0.022803j, +0.683002+0.473351j]]
  wire 3 (     B1): Rot(phi=+0.325283, theta=+1.087839, omega=+1.511560)
      U = [[+0.519466-0.679967j, -0.429101+0.289262j],
       [+0.429101+0.289262j, +0.519466+0.679967j]]
  wire 4 (     B2): Rot(phi=-2.309689, theta=+1.324789, omega=+0.212803)
      U = [[+0.393410+0.683370j, -0.187350+0.585776j],
       [+0.187350+0.585776j, +0.393410-0.683370j]]
```
then **CNOT(control=2, target=4)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [2, 4].

**Block 4** (before CNOT (0, 3)):
```
  wire 0 (a (anc)): Rot(phi=-3.309642, theta=+1.152731, omega=-2.781264)
      U = [[-0.834577+0.080484j, -0.526072+0.142309j],
       [+0.526072+0.142309j, -0.834577-0.080484j]]
  wire 1 (     A1): Rot(phi=-0.441053, theta=+1.986377, omega=+0.017187)
      U = [[+0.533805+0.114856j, -0.815877+0.190275j],
       [+0.815877+0.190275j, +0.533805-0.114856j]]
  wire 2 (     A2): Rot(phi=-0.481848, theta=+1.987561, omega=-1.076572)
      U = [[+0.388124+0.383350j, -0.801313-0.245561j],
       [+0.801313-0.245561j, +0.388124-0.383350j]]
  wire 3 (     B1): Rot(phi=+3.662036, theta=+0.319956, omega=+1.548464)
      U = [[-0.848607-0.504470j, -0.078322-0.138712j],
       [+0.078322-0.138712j, -0.848607+0.504470j]]
  wire 4 (     B2): Rot(phi=+2.166044, theta=+1.214352, omega=+2.742520)
      U = [[-0.634799-0.521059j, -0.547013+0.162187j],
       [+0.547013+0.162187j, -0.634799+0.521059j]]
```
then **CNOT(control=0, target=3)** followed by depolarizing `(1-e2)rho+e2 I/4` on wires [0, 3].

**Block 5** (before measurement):
```
  wire 0 (a (anc)): Rot(phi=+1.195304, theta=+1.047137, omega=+3.008710)
      U = [[-0.438717-0.746695j, -0.308157+0.393717j],
       [+0.308157+0.393717j, -0.438717+0.746695j]]
  wire 1 (     A1): Rot(phi=-3.194684, theta=+2.448653, omega=+2.684158)
      U = [[+0.328576+0.085744j, +0.921420+0.188865j],
       [-0.921420+0.188865j, +0.328576-0.085744j]]
  wire 2 (     A2): Rot(phi=+1.722511, theta=+0.257669, omega=-2.482455)
      U = [[+0.920978+0.367821j, +0.065137-0.110742j],
       [-0.065137-0.110742j, +0.920978-0.367821j]]
  wire 3 (     B1): Rot(phi=-2.558638, theta=+1.630563, omega=+0.581157)
      U = [[+0.376938+0.572758j, -0.000654+0.727918j],
       [+0.000654+0.727918j, +0.376938-0.572758j]]
  wire 4 (     B2): Rot(phi=+2.626628, theta=+0.969038, omega=-2.459053)
      U = [[+0.881795-0.074057j, +0.384755-0.262521j],
       [-0.384755-0.262521j, +0.881795+0.074057j]]
```