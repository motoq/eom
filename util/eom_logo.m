function [gxh, tpv] = eom_logo
% Orbit is an EOM generated Matlab/Octave function that
% plots a 3D orbit trace in ECI or ECF coordinates
%
% Outputs:
%   gxh  Graphics handle to new image
%   tpv  Nx7 matrix of time, pos, vel, in ECI or ECF
%        coordinates, units of DU and DU/TU

tpv = [
  0.0000000000000000e+00  2.9999363043971045e-01 -3.7761041959917152e-01 -9.6440541151487735e-01 9.7588336312604018e-01 7.6874457540201013e-01 2.5640993696526121e-03;
  3.7183429977089427e-01  6.4689649247036973e-01 -8.4358703826080170e-02 -9.1185598088173014e-01 8.7655263688070173e-01 7.9428803378935231e-01 2.7238806371981972e-01;
  7.4366859954178854e-01  9.4584531405009686e-01 2.0385567004559885e-01 -7.7124800039082031e-01 7.2867826732729479e-01 7.4671025955702164e-01 4.6890598936962902e-01;
  1.1155028993126828e+00  1.1896080519601815e+00 4.6668872570703662e-01 -5.7372907372976512e-01 5.8608306493873163e-01 6.6429542831889199e-01 5.8117407832275036e-01;
  1.4873371990835771e+00  1.3851027647752507e+00 6.9735047290341634e-01 -3.4639598655512521e-01 4.7028715859180148e-01 5.7705291690971894e-01 6.3401967002755688e-01;
  1.8591714988544714e+00  1.5427434373058517e+00 8.9684009756279914e-01 -1.0647443648503814e-01 3.8174228152293915e-01 4.9769738810695968e-01 6.5227657511037962e-01;
  2.2310057986253655e+00  1.6716705548504862e+00 1.0688458257501845e+00 1.3638886128982153e-01 3.1481854793488817e-01 4.2931652238521545e-01 6.5180399176643167e-01;
  2.6028400983962601e+00  1.7788329860570671e+00 1.2174077004337280e+00 3.7704857348704029e-01 2.6382742586156388e-01 3.7138474857920833e-01 6.4149353948416699e-01;
  2.9746743981671542e+00  1.8692958679158320e+00 1.3461415600452189e+00 6.1283462555736157e-01 2.2437661275118176e-01 3.2239749717381644e-01 6.2616839961233894e-01;
  3.3465086979380483e+00  1.9467314200357808e+00 1.4580791849608494e+00 8.4241778878550810e-01 1.9331998563299788e-01 2.8079715219650841e-01 6.0844514101423752e-01;
  3.7183429977089428e+00  2.0138254266749263e+00 1.5557093574209424e+00 1.0651988849075831e+00 1.6844517652412214e-01 2.4523858695200301e-01 5.8975505781106441e-01;
  4.0901772974798369e+00  2.0725698432981030e+00 1.6410677995639720e+00 1.2809792009606951e+00 1.4819214950253873e-01 2.1462829047461265e-01 5.7088489272333742e-01;
  4.4620115972507310e+00  2.1244650081175420e+00 1.7158263326155072e+00 1.4897797919355660e+00 1.3144855769684455e-01 1.8809567853784842e-01 5.5226434309141015e-01;
  4.8338458970216251e+00  2.1706584160506228e+00 1.7813679500194457e+00 1.6917405503429956e+00 1.1740996106299589e-01 1.6495079081021702e-01 5.3412175141476959e-01;
  5.2056801967925201e+00  2.2120405899589786e+00 1.8388465162453935e+00 1.8870629016104707e+00 1.0548590422527503e-01 1.4464522431713994e-01 5.1657031227971539e-01;
  5.5775144965634142e+00  2.2493121282239166e+00 1.8892332499213276e+00 2.0759768988795710e+00 9.5236695362240387e-02 1.2674029188394295e-01 4.9965686266517872e-01;
  5.9493487963343084e+00  2.2830312707520486e+00 1.9333527143931237e+00 2.2587222534460980e+00 8.6330437431391099e-02 1.1088226322717068e-01 4.8339004050319939e-01;
  6.3211830961052025e+00  2.3136481453862414e+00 1.9719107323741698e+00 2.4355374696753134e+00 7.8513441240372794e-02 9.6783458650143164e-02 4.6775679983180019e-01;
  6.6930173958760966e+00  2.3415297817664213e+00 2.0055161529228105e+00 2.6066537609448779e+00 7.1589560643363956e-02 8.4207883661653382e-02 4.5273221307273220e-01;
  7.0648516956469916e+00  2.3669786328394111e+00 2.0346979456450525e+00 2.7722918135177337e+00 6.5405550094870132e-02 7.2960297064682211e-02 4.3828532911929141e-01;
  7.4366859954178857e+00  2.3902464653529973e+00 2.0599187313168996e+00 2.9326602538662589e+00 5.9840543038582758e-02 6.2877850897539694e-02 4.2438267629367554e-01;
  7.8085202951887798e+00  2.4115449013231842e+00 2.0815855788907616e+00 3.0879551317766403e+00 5.4798389583029787e-02 5.3823654504525836e-02 4.1099034050073091e-01;
  8.1803545949596739e+00  2.4310535057782032e+00 2.1000586905888676e+00 3.2383600012230254e+00 5.0202005523081482e-02 4.5681782648468179e-02 3.9807517297500000e-01;
  8.5521888947305680e+00  2.4489260545418889e+00 2.1156584428164007e+00 3.3840463428671850e+00 4.5989155055864796e-02 3.8353373309138017e-02 3.8560546313982402e-01;
  8.9240231945014621e+00  2.4652954365244906e+00 2.1286711369563198e+00 3.5251741706114710e+00 4.2109268386832828e-02 3.1753553311840196e-02 3.7355128231568163e-01;
  9.2958574942723562e+00  2.4802775204440408e+00 2.1393537299613121e+00 3.6618927253867657e+00 3.8521015304856021e-02 2.5808997485421963e-02 3.6188462579047642e-01;
  9.6676917940432503e+00  2.4939742282826263e+00 2.1479377520521541e+00 3.7943411971607337e+00 3.5190437210415274e-02 2.0455976339489444e-02 3.5057943292398824e-01;
  1.0039526093814146e+01  2.5064759953829086e+00 2.1546325719628223e+00 3.9226494398345455e+00 3.2089496061701833e-02 1.5638783310993151e-02 3.3961153532579236e-01;
  1.0411360393585040e+01  2.5178637521352396e+00 2.1596281348415314e+00 4.0469386585874512e+00 2.9194937673184613e-02 1.1308459137186104e-02 3.2895856459041434e-01;
  1.0783194693355934e+01  2.5282105294657211e+00 2.1630972710845064e+00 4.1673220585910586e+00 2.6487394248228616e-02 7.4217505123687745e-03 3.1859983934999592e-01;
  1.1155028993126828e+01  2.5375827662406634e+00 2.1651976538556892e+00 4.2839054498724893e+00 2.3950670575359878e-02 3.9402547753968299e-03 3.0851624394649702e-01;
  1.1526863292897723e+01  2.5460413787928657e+00 2.1660734672373811e+00 4.3967878067280042e+00 2.1571172388278135e-02 8.2971331016339165e-04 2.9869010625858589e-01;
  1.1898697592668617e+01  2.5536426393502238e+00 2.1658568346805014e+00 4.5060617822809430e+00 1.9337445619877075e-02 -1.9405754046906178e-03 2.8910507916547579e-01;
  1.2270531892439511e+01  2.5604388999834686e+00 2.1646690478390651e+00 4.6118141800543935e+00 1.7239802789598527e-02 -4.3982458862580759e-03 2.7974602817698491e-01;
  1.2642366192210405e+01  2.5664791909400417e+00 2.1626216283316322e+00 4.7141263851206867e+00 1.5270018325239510e-02 -6.5682768650582417e-03 2.7059892651623557e-01;
  1.3014200491981299e+01  2.5718097162760425e+00 2.1598172489984084e+00 4.8130747577116146e+00 1.3421078775076331e-02 -8.4733379162148324e-03 2.6165075815892280e-01;
  1.3386034791752193e+01  2.5764742650814281e+00 2.1563505364710216e+00 4.9087309922630853e+00 1.1686976995174287e-02 -1.0134082919798258e-02 2.5288942884978127e-01;
  1.3757869091523087e+01  2.5805145529984559e+00 2.1523087730581367e+00 5.0011624448172514e+00 1.0062541771321662e-02 -1.1569399525223321e-02 2.4430368482788417e-01;
  1.4129703391293983e+01  2.5839705059104667e+00 2.1477725128785190e+00 5.0904324315719451e+00 8.5432961503006458e-03 -1.2796622034680245e-02 2.3588303882922199e-01;
  1.4501537691064877e+01  2.5868804954468816e+00 2.1428161246861182e+00 5.1766005011893625e+00 7.1253391525143662e-03 -1.3831713730134769e-02 2.2761770285273594e-01;
  1.4873371990835771e+01  2.5892815341800284e+00 2.1375082718036187e+00 5.2597226832780066e+00 5.8052466206597320e-03 -1.4689423563729020e-02 2.1949852714494572e-01;
  1.5245206290606665e+01  2.5912094369735734e+00 2.1319123379220790e+00 5.3398517152586820e+00 4.5799878032556579e-03 -1.5383421248830565e-02 2.1151694485921332e-01;
  1.5617040590377560e+01  2.5926989538062033e+00 2.1260868061595213e+00 5.4170372496263042e+00 3.4468549339102816e-03 -1.5926414080287794e-02 2.0366492186560939e-01;
  1.5988874890148454e+01  2.5937838784761134e+00 2.1200855976437798e+00 5.4913260434298516e+00 2.4034035893989981e-03 -1.6330248240378063e-02 1.9593491121808199e-01;
  1.6360709189919348e+01  2.5944971368476399e+00 2.1139583749493145e+00 5.5627621316161822e+00 1.4474020237369965e-03 -1.6605996883043369e-02 1.8831981182165494e-01;
  1.6732543489690244e+01  2.5948708576950739e+00 2.1077508149372743e+00 5.6313869857201162e+00 5.7678800549959899e-04 -1.6764036911006492e-02 1.8081293088036587e-01;
  1.7104377789461136e+01  2.5949364287007128e+00 2.1015048548972004e+00 5.6972396592346586e+00 -2.1036804998103475e-04 -1.6814116051015005e-02 1.7340794974434925e-01;
  1.7476212089232032e+01  2.5947245397579404e+00 2.0952589153381536e+00 5.7603569208598477e+00 -9.1589464918535543e-04 -1.6765411578198667e-02 1.6609889281065757e-01;
  1.7848046389002924e+01  2.5942652153913337e+00 2.0890481023169896e+00 5.8207733767068666e+00 -1.5415463134372275e-03 -1.6626581830695601e-02 1.5888009916634385e-01;
  1.8219880688773820e+01  2.5935878378262940e+00 2.0829043918005077e+00 5.8785215824240620e+00 -2.0890268805318734e-03 -1.6405811481823129e-02 1.5174619669363271e-01;
  1.8591714988544712e+01  2.5927211620096937e+00 2.0768567982245871e+00 5.9336321461127062e+00 -2.5600095853876106e-03 -1.6110851392492697e-02 1.4469207838560755e-01;
  1.8963549288315608e+01  2.5916933236876698e+00 2.0709315291316752e+00 5.9861338228117535e+00 -2.9561543971245098e-03 -1.5749053745824484e-02 1.3771288064672579e-01;
  1.9335383588086501e+01  2.5905318414831857e+00 2.0651521275274471e+00 6.0360536012513260e+00 -3.2791230136101539e-03 -1.5327403064832608e-02 1.3080396337576647e-01;
  1.9707217887857396e+01  2.5892636137813012e+00 2.0595396033871811e+00 6.0834167835032957e+00 -3.5305918516427156e-03 -1.4852543629038725e-02 1.2396089164971670e-01;
  2.0079052187628292e+01  2.5879149111103246e+00 2.0541125555692390e+00 6.1282470580934474e+00 -3.7122633189587582e-03 -1.4330803734209108e-02 1.1717941884574587e-01;
  2.0450886487399185e+01  2.5865113646144815e+00 2.0488872852325479e+00 6.1705665670823162e+00 -3.8258756112472805e-03 -1.3768217178784063e-02 1.1045547105502042e-01;
  2.0822720787170081e+01  2.5850779511244912e+00 2.0438779017291289e+00 6.2103959675701059e+00 -3.8732112415505822e-03 -1.3170542309122367e-02 1.0378513265692429e-01;
  2.1194555086940973e+01  2.5836389752641673e+00 2.0390964218237952e+00 6.2477544880348432e+00 -3.8561044796922104e-03 -1.2543278911921085e-02 9.7164632935321960e-02;
  2.1566389386711869e+01  2.5822180489699829e+00 2.0345528629927614e+00 6.2826599798710330e+00 -3.7764478545320437e-03 -1.1891683204807602e-02 9.0590333630162514e-02;
  2.1938223686482761e+01  2.5808380687457202e+00 2.0302553314682497e+00 6.3151289644583652e+00 -3.6361978511105725e-03 -1.1220781144163812e-02 8.4058717327986948e-02;
  2.2310057986253657e+01  2.5795211909313069e+00 2.0262101056164368e+00 6.3451766760559520e+00 -3.4373799174039800e-03 -1.0535380241818697e-02 7.7566376604035453e-02;
  2.2681892286024549e+01  2.5782888052216801e+00 2.0224217151735355e+00 6.3728171007865315e+00 -3.1820928809026330e-03 -9.8400800586887251e-03 7.1110003836659302e-02;
  2.3053726585795445e+01  2.5771615066396967e+00 2.0188930168023416e+00 6.3980630119470607e+00 -2.8725128631241696e-03 -9.1392815231399853e-03 6.4686381621848796e-02;
  2.3425560885566338e+01  2.5761590661345206e+00 2.0156252663820662e+00 6.4209260018563290e+00 -2.5108967700499260e-03 -8.4371952043010637e-03 5.8292373721922634e-02;
  2.3797395185337233e+01  2.5753003999506405e+00 2.0126181883974081e+00 6.4414165104268859e+00 -2.0995854280822903e-03 -7.7378486553894271e-03 5.1924916487891118e-02;
  2.4169229485108129e+01  2.5746035378859489e+00 2.0098700427553431e+00 6.4595438506268392e+00 -1.6410064281585910e-03 -7.0450929289692585e-03 4.5581010699808820e-02;
  2.4541063784879022e+01  2.5740855905384281e+00 2.0073776893180577e+00 6.4753162309773122e+00 -1.1376767349781261e-03 -6.3626083546287687e-03 3.9257713773599948e-02;
  2.4912898084649918e+01  2.5737627156178360e+00 2.0051366504124473e+00 6.4887407752127775e+00 -5.9220511367985378e-04 -5.6939096596618794e-03 3.2952132286450622e-02;
  2.5284732384420810e+01  2.5736500833813638e+00 2.0031411715464746e+00 6.4998235392143426e+00 -7.2944226498793493e-06 -5.0423505046964271e-03 2.6661414775986242e-02;
  2.5656566684191706e+01  2.5737618412348522e+00 2.0013842805373154e+00 6.5085695253098130e+00 6.1425618166366012e-04 -4.4111274987255296e-03 2.0382744771084850e-02;
  2.6028400983962598e+01  2.5741110775230824e+00 1.9998578452350693e+00 6.5149826940192890e+00 1.2695490841135369e-03 -3.8032837514504216e-03 1.4113334014399605e-02;
  2.6400235283733494e+01  2.5747097845189266e+00 1.9985526300026448e+00 6.5190659733103038e+00 1.9555860519240372e-03 -3.2217120151921282e-03 7.8504158384885368e-03;
  2.6772069583504386e+01  2.5755688206031864e+00 1.9974583510967754e+00 6.5208212654128879e+00 2.6692660357174429e-03 -2.6691574637052024e-03 1.5912386588953110e-03;
  2.7143903883275282e+01  2.5766978716133067e+00 1.9965637310759878e+00 6.5202494512313613e+00 3.4073828888271443e-03 -2.1482201509927090e-03 -4.6669404513767542e-03;
  2.7515738183046174e+01  2.5781054113214203e+00 1.9958565523488403e+00 6.5173503923765841e+00 4.1666229627389190e-03 -1.6613571895764268e-03 -1.0926862140784165e-02;
  2.7887572482817070e+01  2.5797986609899985e+00 1.9953237099568113e+00 6.5121229308296948e+00 4.9435625353896244e-03 -1.2108846845943167e-03 -1.7191271145924121e-02;
  2.8259406782587966e+01  2.5817835479318560e+00 1.9949512636800404e+00 6.5045648862354346e+00 5.7346650273862529e-03 -7.9897945749220862e-04 -2.3462922680209500e-02;
  2.8631241082358859e+01  2.5840646629891579e+00 1.9947244895351564e+00 6.4946730508105492e+00 6.5362779588664438e-03 -4.2768059095876109e-04 -2.9744588878993367e-02;
  2.9003075382129754e+01  2.5866452168249729e+00 1.9946279307272858e+00 6.4824431818397379e+00 7.3446295967909006e-03 -9.8890825070049456e-05 -3.6039065336109379e-02;
  2.9374909681900647e+01  2.5895269949015143e+00 1.9946454481070877e+00 6.4678699917186098e+00 8.1558252387737371e-03 1.8562216666738449e-04 -4.2349177767652964e-02;
  2.9746743981671543e+01  2.5927103110011176e+00 1.9947602701700118e+00 6.4509471354893915e+00 8.9658430751497209e-03 4.2422459349643511e-04 -4.8677788840027533e-02;
  3.0118578281442435e+01  2.5961939591198995e+00 1.9949550426288976e+00 6.4316671958012979e+00 9.7705295657075988e-03 6.1541542837692705e-04 -5.5027805200813962e-02;
  3.0490412581213331e+01  2.5999751635411221e+00 1.9952118775786336e+00 6.4100216652126081e+00 1.0565594261325393e-02 7.5782532626950919e-04 -6.1402184752970310e-02;
  3.0862246880984223e+01  2.6040495268686028e+00 1.9955124022605002e+00 6.3860009257361066e+00 1.1346603993490363e-02 8.5021555178225430e-04 -6.7803944215178866e-02;
  3.1234081180755119e+01  2.6084109757685754e+00 1.9958378074255920e+00 6.3595942255131268e+00 1.2108976346236677e-02 8.9147688724741192e-04 -7.4236167013942686e-02;
  3.1605915480526011e+01  2.6130517041366583e+00 1.9961688952832308e+00 6.3307896524838494e+00 1.2847972315205547e-02 8.8062849076271416e-04 -8.0702011556255809e-02;
  3.1977749780296907e+01  2.6179621133668101e+00 1.9964861270117149e+00 6.2995741049026952e+00 1.3558688047146460e-02 8.1681667159114819e-04 -8.7204719935441932e-02;
  3.2349584080067800e+01  2.6231307493606941e+00 1.9967696697912503e+00 6.2659332585273120e+00 1.4236045539948053e-02 6.9931354745298346e-04 -9.3747627127092142e-02;
  3.2721418379838696e+01  2.6285442358637483e+00 1.9969994433113563e+00 6.2298515302875588e+00 1.4874782167956667e-02 5.2751554461712917e-04 -1.0033417073697649e-01;
  3.3093252679609591e+01  2.6341872036661482e+00 1.9971551656825972e+00 6.1913120382169424e+00 1.5469438879532262e-02 3.0094169707659712e-04 -1.0696790136847420e-01;
  3.3465086979380487e+01  2.6400422151399425e+00 1.9972163986734994e+00 6.1502965574024655e+00 1.6014346893103160e-02 1.9231695428818511e-05 -1.1365249368353374e-01;
  3.3836921279151376e+01  2.6460896835214163e+00 1.9971625921647282e+00 6.1067854716803138e+00 1.6503612693874240e-02 -3.1785637095730523e-04 -1.2039175823844399e-01;
  3.4208755578922272e+01  2.6523077862635072e+00 1.9969731276979308e+00 6.0607577207727230e+00 1.6931101105274984e-02 -7.1044864346727808e-04 -1.2718965418403011e-01;
  3.4580589878693168e+01  2.6586723716968446e+00 1.9966273609651499e+00 6.0121907425266343e+00 1.7290416176377048e-02 -1.1585585708975788e-03 -1.3405030292923936e-01;
  3.4952424178464064e+01  2.6651568581313234e+00 1.9961046630594734e+00 5.9610604098758646e+00 1.7574879588083157e-02 -1.6620905214057343e-03 -1.4097800287767404e-01;
  3.5324258478234952e+01  2.6717321244161627e+00 1.9953844602655542e+00 5.9073409621055628e+00 1.7777506235746549e-02 -2.2208440334335732e-03 -1.4797724535860565e-01;
  3.5696092778005848e+01  2.6783663908348356e+00 1.9944462721326905e+00 5.8510049299502516e+00 1.7890976592722312e-02 -2.8345188280914017e-03 -1.5505273188743349e-01;
  3.6067927077776744e+01  2.6850250890554919e+00 1.9932697475205163e+00 5.7920230540036846e+00 1.7907605396682895e-02 -3.5027207284919165e-03 -1.6220939290571598e-01;
  3.6439761377547640e+01  2.6916707196729770e+00 1.9918346982485924e+00 5.7303641958597558e+00 1.7819306126232344e-02 -4.2249686595273546e-03 -1.6945240816793194e-01;
  3.6811595677318536e+01  2.6982626956644369e+00 1.9901211299112436e+00 5.6659952413382912e+00 1.7617550647240878e-02 -5.0007029358492975e-03 -1.7678722896112639e-01;
  3.7183429977089425e+01  2.7047571698321655e+00 1.9881092693321976e+00 5.5988809950762484e+00 1.7293323303290913e-02 -5.8292950878276887e-03 -1.8421960236491905e-01;
  3.7555264276860321e+01  2.7111068440128601e+00 1.9857795880341216e+00 5.5289840656834208e+00 1.6837068599117027e-02 -6.7100595268638633e-03 -1.9175559778299536e-01;
  3.7927098576631217e+01  2.7172607574889960e+00 1.9831128209746849e+00 5.4562647405711200e+00 1.6238631475412733e-02 -7.6422674150816100e-03 -1.9940163600337479e-01;
  3.8298932876402112e+01  2.7231640516365809e+00 1.9800899796457323e+00 5.3806808494612772e+00 1.5487188992153266e-02 -8.6251631831827571e-03 -2.0716452107351999e-01;
  3.8670767176173001e+01  2.7287577073581235e+00 1.9766923584584026e+00 5.3021876154712508e+00 1.4571172018699624e-02 -9.6579842379452847e-03 -2.1505147530762839e-01;
  3.9042601475943897e+01  2.7339782512873416e+00 1.9729015331043693e+00 5.2207374925455809e+00 1.3478175263311340e-02 -1.0739984522715133e-02 -2.2307017777703719e-01;
  3.9414435775714793e+01  2.7387574260672807e+00 1.9686993493190557e+00 5.1362799878688659e+00 1.2194853651269832e-02 -1.1870462746631833e-02 -2.3122880666997667e-01;
  3.9786270075485689e+01  2.7430218191989200e+00 1.9640679001254795e+00 5.0487614677435477e+00 1.0706802665016561e-02 -1.3048796289978737e-02 -2.3953608594307249e-01;
  4.0158104375256585e+01  2.7466924439768969e+00 1.9589894892278803e+00 4.9581249452527834e+00 8.9984197734975756e-03 -1.4274482034925838e-02 -2.4800133672221764e-01;
  4.0529938675027473e+01  2.7496842648578377e+00 1.9534465776981267e+00 4.8643098478523310e+00 7.0527434775209275e-03 -1.5547185677911982e-02 -2.5663453394216301e-01;
  4.0901772974798369e+01  2.7519056581764350e+00 1.9474217104572866e+00 4.7672517628486393e+00 4.8512657532342029e-03 -1.6866801471323228e-02 -2.6544636873808164e-01;
  4.1273607274569265e+01  2.7532577973969996e+00 1.9408974182309544e+00 4.6668821585272013e+00 2.3737127471245676e-03 -1.8233524844080581e-02 -2.7444831711194334e-01;
  4.1645441574340161e+01  2.7536339499590912e+00 1.9338560896410948e+00 4.5631280785025519e+00 -4.0221258811618085e-04 -1.9647940997946423e-02 -2.8365271538176734e-01;
  4.2017275874111050e+01  2.7529186701724009e+00 1.9262798067890581e+00 4.4559118066808310e+00 -3.5011336949569262e-03 -2.1111133416062328e-02 -2.9307284286739677e-01;
  4.2389110173881946e+01  2.7509868693885617e+00 1.9181501360306201e+00 4.3451505000767856e+00 -6.9503590881446380e-03 -2.2624817316450212e-02 -3.0272301214917618e-01;
  4.2760944473652842e+01  2.7477027406731027e+00 1.9094478635132905e+00 4.2307557866383982e+00 -1.0780308519260741e-02 -2.4191504524053914e-02 -3.1261866702070684e-01;
  4.3132778773423738e+01  2.7429185101894733e+00 1.9001526623048137e+00 4.1126333252487717e+00 -1.5025021686718664e-02 -2.5814708142145197e-02 -3.2277648789043634e-01;
  4.3504613073194626e+01  2.7364729812099511e+00 1.8902426743724356e+00 3.9906823252701336e+00 -1.9722768485491940e-02 -2.7499197947359191e-02 -3.3321450378752826e-01;
  4.3876447372965522e+01  2.7281898286934747e+00 1.8796939860093784e+00 3.8647950234715682e+00 -2.4916784912471211e-02 -2.9251320851682465e-02 -3.4395220917146030e-01;
  4.4248281672736418e+01  2.7178755922144258e+00 1.8684799691517149e+00 3.7348561171115446e+00 -3.0656165454861203e-02 -3.1079405409204834e-02 -3.5501068223977206e-01;
  4.4620115972507314e+01  2.7053173019925212e+00 1.8565704528572498e+00 3.6007421535797204e+00 -3.6996951650090287e-02 -3.2994275683600478e-02 -3.6641269907764024e-01;
  4.4991950272278210e+01  2.6902796559297251e+00 1.8439306782669898e+00 3.4623208797367391e+00 -4.4003468295165334e-02 -3.5009908542100121e-02 -3.7818283433106054e-01;
  4.5363784572049099e+01  2.6725016436252993e+00 1.8305199755616681e+00 3.3194505585379654e+00 -5.1749974594479962e-02 -3.7144280641741814e-02 -3.9034753337793188e-01;
  4.5735618871819995e+01  2.6516924845479553e+00 1.8162900812103928e+00 3.1719792676499909e+00 -6.0322718912869056e-02 -3.9420468564670373e-02 -4.0293513204485848e-01;
  4.6107453171590890e+01  2.6275267094382353e+00 1.8011829859333395e+00 3.0197442060878377e+00 -6.9822514944100342e-02 -4.1868090055378962e-02 -4.1597578587565542e-01;
  4.6479287471361786e+01  2.5996381631610519e+00 1.7851281649341140e+00 2.8625710528473034e+00 -8.0367997162115803e-02 -4.4525209633198844e-02 -4.2950124869408357e-01;
  4.6851121771132675e+01  2.5676126387749965e+00 1.7680389871478359e+00 2.7002734500165273e+00 -9.2099768918083066e-02 -4.7440883420960679e-02 -4.4354440455542404e-01;
  4.7222956070903571e+01  2.5309787596583893e+00 1.7498080219923688e+00 2.5326527283781157e+00 -1.0518573399371228e-01 -5.0678594314349720e-02 -4.5813839943143175e-01;
  4.7594790370674467e+01  2.4891965993292211e+00 1.7303008489231029e+00 2.3594980668056915e+00 -1.1982801115924162e-01 -5.4320943024418047e-02 -4.7331512414528210e-01;
  4.7966624670445363e+01  2.4416433531717336e+00 1.7093478091048848e+00 2.1805873960183386e+00 -1.3627198453371689e-01 -5.8476134572237957e-02 -4.8910264195739794e-01;
  4.8338458970216259e+01  2.3875951328626512e+00 1.6867328915995567e+00 1.9956895538649093e+00 -1.5481825846158001e-01 -6.3287068465311486e-02 -5.0552088600557510e-01;
  4.8710293269987147e+01  2.3262036154031724e+00 1.6621785736319517e+00 1.8045685284591251e+00 -1.7583858746242303e-01 -6.8944261484516142e-02 -5.2257448818160879e-01;
  4.9082127569758043e+01  2.2564658076422557e+00 1.6353248630375365e+00 1.6069911863283832e+00 -1.9979726440150797e-01 -7.5704499816139845e-02 -5.4024078281675236e-01;
  4.9453961869528939e+01  2.1771845399332483e+00 1.6056999022043084e+00 1.4027408574553228e+00 -2.2727998414396219e-01 -8.3918189412011096e-02 -5.5844955197868307e-01;
  4.9825796169299835e+01  2.0869164406050591e+00 1.5726780924942978e+00 1.1916408803665470e+00 -2.5903280051759547e-01 -9.4070105913304766e-02 -5.7704835109182573e-01;
  5.0197630469070724e+01  1.9839030818920091e+00 1.5354194719883709e+00 9.7359535810201414e-01 -2.9601418304370802e-01 -1.0684103379973374e-01 -5.9574209129572597e-01;
  5.0569464768841620e+01  1.8659799518168312e+00 1.4927805397203771e+00 7.4866023803453041e-01 -3.3946239057989025e-01 -1.2320215143861400e-01 -6.1398556417019778e-01;
  5.0941299068612516e+01  1.7304577916773665e+00 1.4431812004691638e+00 5.1716899840175168e-01 -3.9097539771032835e-01 -1.4456028583677638e-01 -6.3078792145260776e-01;
  5.1313133368383411e+01  1.5739745859769609e+00 1.3844044762948544e+00 2.7995888443738542e-01 -4.5258239412797546e-01 -1.7297874588877188e-01 -6.4434924946998640e-01;
  5.1684967668154307e+01  1.3923330960785920e+00 1.3132964341883522e+00 3.8785749948303083e-02 -5.2672788564628847e-01 -2.1149476786913218e-01 -6.5137502522498802e-01;
  5.2056801967925196e+01  1.1803947058372797e+00 1.2253339076559917e+00 -2.0290555340664734e-01 -6.1591679910851815e-01 -2.6449429462513613e-01 -6.4579085332914632e-01;
  5.2428636267696092e+01  9.3226810923382331e-01 1.1140836244603152e+00 -4.3856258567972617e-01 -7.2129616039601230e-01 -3.3781597533543006e-01 -6.1648554928139521e-01;
  5.2800470567466988e+01  6.4247454460694997e-01 9.7085335899444802e-01 -6.5608936324329226e-01 -8.3841491571060256e-01 -4.3723255989155868e-01 -5.4427992435451178e-01;
  5.3172304867237884e+01  3.0961134315982097e-01 7.8574329502024443e-01 -8.3464090167806648e-01 -9.4778794822003076e-01 -5.6159157571531937e-01 -4.0217937554652317e-01;
  5.3544139167008773e+01  -5.6045657691822189e-02 5.5287392118807366e-01 -9.4418582432033094e-01 -1.0058404465479187e+00 -6.8755848284086796e-01 -1.7327222796453179e-01
];

n = size(tpv,1);

acc = tpv(:,2:4);
for ii = 1:n
  rmag = norm(acc(ii,:));
  acc(ii,:) = -acc(ii,:)/(rmag*rmag*rmag);
end

gxh = figure; hold on;
%scatter3(0, 0, 0, 'b');

plot3(acc(:,1), acc(:,2), acc(:,3), '.r');%, 'MarkerSize', 8);
%scatter3(acc(1,1), acc(1,2), acc(1,3), 'g');
%scatter3(acc(n,1), acc(n,2), acc(n,3), 'r');

plot3(tpv(:,5), tpv(:,6), tpv(:,7), '.m');%, 'MarkerSize', 8);
%scatter3(tpv(1,5), tpv(1,6), tpv(1,7), 'g');
%scatter3(tpv(n,5), tpv(n,6), tpv(n,7), 'r');

plot3(tpv(:,2), tpv(:,3), tpv(:,4), '.b');%, 'MarkerSize', 8);
%scatter3(tpv(1,2), tpv(1,3), tpv(1,4), 'g');
%scatter3(tpv(n,2), tpv(n,3), tpv(n,4), 'r');
axis off;
axis equal;
%view([40 0])
view([220 0]);

end
