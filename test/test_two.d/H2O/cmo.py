import numpy
from util.full import matrix

cmo = matrix((24, 24))
cmo[:11, :11] = numpy.array([
    [ 1.00037005644774, 0.00213928142197, -0.00139645162733,  0.00250320614892, -0.00168365551156, 0.00000417179213, -0.00010568751627, -0.00027311036476,  0.00056204660778, 0.00052477299313,  0.00045051613113],
    [-0.00764995676695, 0.87372007927195, -0.11446369166089,  0.11081294767437, -0.06188049075196, 0.00122390477237, -0.00288482651200,  0.31909593297047, -0.14943668512171,-0.03635107656087, -0.02008905901125],
    [ 0.00218188411380,-0.29318297900759, -0.17682717408706,  0.79267177820451,  0.00181832874368, 0.01768753372185, -0.00458766924852,  0.35155848740093, -0.13271279540296,-0.03165237343243,  0.00704103155682],
    [ 0.05189274589511,-0.14332375175701, -0.90584515407463, -0.28547304713507, -0.21170912054705,-0.00753222886420,  0.00920238615519,  0.10575589422514,  0.76987725358803,-0.01621298099706, -0.01563551317394],
    [ 0.06833791423126, 0.50854763250820, -0.61069110991343,  0.47742207558481, -0.21760029506920,-0.00037492684871,  0.10481322637675, -1.36310352174187,  1.24086084449933,-0.28378020498683, -0.05694956930942],
    [ 0.00307036493807, 0.22313096119853,  0.58194482317608, -1.10598495943482,  1.76491922320823, 0.01456983583454,  0.00741605994232, -0.88418613348653,  0.33025487687161, 0.07966186404158,  0.26045387900112],
    [ 0.20385065017208, 1.04814507440774, -1.09392840725873, -0.21723508898129,  0.84988495567017,-0.12113033846152, -0.00934519422054, -0.57077514515853,  0.49008704403139, 0.33901487411969, -0.53541517404592],
    [-0.73049680997176,-3.16089341650320,  4.80229329699141,  0.12075736306732,  0.91054393821628,-0.11542286282910,  0.10877087324223, -1.37094247910108,  0.30453649173373, 0.37030953438007, -0.02406337819401],
    [ 0.24567681999547, 1.31174602774700,  0.07480159732451,  0.95885257738136, -0.24905823475546,-0.05207781435411,  0.20872718389277, -0.79551074266912,  0.23277611879307, 0.74469921902456,  0.53405215283822],
    [-0.08304755903324,-0.30102526839298,  0.69763628042051,  0.01161994402264,  0.51175312137587, 1.07557220332225,  0.18150224264537, -0.48808704111997,  0.20926860411269, 0.34351066042799, -0.30977269672795],
    [ 0.16542053277834, 0.43728192076403, -2.47184521157967, -0.55149357116402, -0.65786651435936,-0.14498630937364,  1.15464174017908,  2.02553837169866, -0.82422195869892,-0.61616382031630, -0.48136846319370]
    ])

cmo[11:15, 11:15] = numpy.array([
    [ 0.92041559956108, 0.07241994018821, 0.01799623382026, 0.03047987898719],
    [-1.40987905655648, 1.66245717885341,-0.00243195680743, 0.02286179295284],
    [-0.01140426393900, 0.64793492608957,-0.16187699280932,-0.76395012286075],
    [ 0.02640107627798,-0.31339796099761,-1.03853627144568, 0.38715704931552]
    ])

cmo[15:22, 15:22] = numpy.array([
    [ 0.71107127316501,-0.10120394303449, 0.02623656642151, 0.55299314427539, -0.18452319045432,-0.02341191289942,-0.03243341057063],
    [-0.41903420742831,-0.47748790342782,-0.02012495853895, 0.05011639045556,  1.39118782976748,-0.01977454415704,-0.01572940892190],
    [-0.39661654197134,-0.25149402388002, 0.11029302744940, 1.60368717121215, -1.51258385409427, 0.06999356017923, 0.14852461859716],
    [-1.07852596404538, 2.24437098774909, 0.03997647831412,-0.63075245833921, -0.48585979264347, 0.30586093347818, 0.17872482219478],
    [ 0.03359009460491,-0.91572119653225, 0.02803938305135, 0.65513796399866, -0.23664365611728, 0.45883057276844,-0.69301562063371],
    [ 1.19613528805294,-0.39762263698350,-0.13142659766418,-0.46262097822477, -0.10003272808805, 0.71647466571577, 0.54825726982745],
    [ 0.63634507688857, 0.83117009927124, 1.31829718409079,-1.76715139017784,  0.60175120616249, 0.59389397795861, 0.48771483717728]
    ])

cmo[22:24, 22:24] = numpy.array([
    [-0.13336003346596 -0.68435603395917],
    [ 1.06074245612680 -0.36200250458266]
    ])

if __name__ == "__main__":
    print cmo
