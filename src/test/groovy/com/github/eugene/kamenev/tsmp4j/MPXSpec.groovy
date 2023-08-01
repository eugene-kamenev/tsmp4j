package com.github.eugene.kamenev.tsmp4j

class MPXSpec extends BaseSpec {
  def 'test moving stats method'() {
    given:
    var meanx = new double[]{56.97833333333333, 63.88666666666666, 65.14166666666667, 66.98666666666668, 68.92833333333334,
      67.07000000000001, 62.610000000000014, 55.49666666666669, 51.97166666666669, 50.76833333333335, 48.416666666666686,
      47.23333333333335, 46.95000000000001, 48.238333333333344, 47.623333333333335, 43.80666666666665, 38.54999999999999,
      35.31166666666667, 33.916666666666686, 31.30666666666669, 30.18833333333335, 32.00833333333336, 35.246666666666705,
      37.67833333333336, 38.964999999999996, 43.554999999999986, 47.28666666666667, 46.54, 47.47666666666665,
      51.72333333333332, 52.291666666666664, 49.10499999999998, 47.04499999999996, 44.2483333333333, 39.19999999999997,
      29.306666666666654, 21.85833333333331, 19.79833333333333, 16.660000000000007, 15.941666666666682, 16.958333333333332,
      20.49666666666667, 23.42500000000003, 23.00333333333337, 24.860000000000014, 26.518333333333317, 27.37833333333333,
      29.634999999999952, 38.554999999999914, 48.34999999999991, 55.056666666666615, 62.63166666666666, 72.5,
      80.98333333333335, 81.76500000000003, 82.06666666666668, 83.31333333333335, 81.6583333333333, 76.90166666666664,
      72.91000000000001, 70.1833333333334, 65.42500000000003, 60.30000000000003, 56.98166666666672, 54.22666666666669,
      49.63333333333336, 47.32166666666664, 47.166666666666664, 50.086666666666666, 51.685000000000024, 52.97833333333339,
      54.7166666666667, 52.081666666666706, 48.49333333333334, 44.98333333333335, 44.39499999999998, 43.714999999999996,
      41.48000000000002, 43.53333333333338, 47.62500000000008, 48.433333333333394, 47.84833333333336, 46.823333333333345,
      47.24833333333337, 44.24666666666675, 39.691666666666755, 33.60666666666672, 31.19333333333346, 28.82500000000012,
      26.490000000000084, 25.428333333333438, 27.071666666666715, 28.783333333333456, 27.95500000000008, 29.285000000000007,
      33.61333333333338, 37.22166666666666, 36.64000000000002, 38.559999999999945, 40.916666666666664, 41.77666666666679,
      39.26333333333347, 38.36166666666683, 40.30166666666688, 39.968333333333554, 35.981666666666875, 33.206666666666784,
      30.625, 27.283333333333303, 27.221666666666653, 29.171666666666624, 36.87833333333325, 45.43166666666669,
      52.5916666666667, 59.88833333333332, 62.566666666666606, 66.02333333333338, 63.04333333333337, 57.964999999999996,
      55.30000000000003, 54.41166666666671, 53.52333333333339, 49.24166666666664}
    var stdx = new double[] {0.03643288531388835, 0.0460826491969164, 0.05044988207594219, 0.05654416289195558,
      0.06995123341670403, 0.045554838854224955, 0.0382909684546397, 0.03319218068464188, 0.032294737288699915,
      0.036011815729799966, 0.04886237935555518, 0.0477463496512956, 0.04770538487014686, 0.05720972808418981,
      0.052465497605753045, 0.07095266091722408, 0.05587748528648964, 0.04629343463103532, 0.05378047082628052,
      0.060647516315919005, 0.06802811073303973, 0.0455013172396381, 0.046314240717158166, 0.0577965278531788,
      0.053560262096294886, 0.05510568262243995, 0.06909084334234278, 0.06854574747208356, 0.07100379510753035,
      0.06219885401292658, 0.06668959948196139, 0.0477730824821783, 0.047354156908567244, 0.035974246836774186,
      0.02797384142628936, 0.027451736315658493, 0.02827554920874882, 0.031018019040144864, 0.04408450366162548,
      0.047868402015141774, 0.04492553913821344, 0.05524589261654207, 0.34573008784547404, 0.259598145048428,
      0.09423451031750883, 0.07850012220238016, 0.07857094054451293, 0.061908201520319985, 0.02360785974966055,
      0.019724504318366868, 0.0190486068939605, 0.020052300907551814, 0.025949580370004045, 0.06642702673008384,
      0.07244833397361812, 0.0735543238147497, 0.09418570610146569, 0.05889423644901258, 0.042086197105371245,
      0.04860393754367375, 0.04887543889931355, 0.0463258874515598, 0.06219037760016283, 0.06259740946446586,
      0.051468614882874686, 0.052119544753828725, 0.08347743212109471, 0.08631590800471516, 0.04387709347632502,
      0.04185554634862731, 0.044895167175035146, 0.05965315881910853, 0.039228107207148276, 0.030867391976453978,
      0.04146846129856282, 0.044835199414318946, 0.04692326465330002, 0.04617603456110419, 0.04956921118892928,
      0.06243514138244356, 0.06050328305973241, 0.06474724880240418, 0.06013248310605006, 0.06690364526629235,
      0.04330380353583133, 0.04154376484021352, 0.03582729402321518, 0.047629355827931465, 0.05542379413824536,
      0.06521785054319207, 0.06351275664897432, 0.05205921272042949, 0.06704314546945747, 0.0783306224281164,
      0.07068104246939184, 0.04557399790396652, 0.052759475961726136, 0.05207799357844761, 0.06557564561904113,
      0.06901110205339461, 0.07658819565493995, 0.09265030518760453, 0.10211011644638356, 0.09837999257182915,
      0.09203387788055925, 0.05718184900839983, 0.05001808897826772, 0.042592591953557656, 0.04116612741432197,
      0.041639628361947574, 0.033645701727296064, 0.022303682590818247, 0.018808638233379787, 0.021730130708182055,
      0.03827389348674806, 0.04864744729540487, 0.08405892417029395, 0.05543162874745679, 0.050997092704192976,
      0.04658153503101733, 0.049035677429447096, 0.05216280877327378, 0.07467298217497143}

    when:
    var stats = MPX.muinvn(timeSeries, 6)

    then:
    equals(meanx, stats[0])
    equals(stdx, stats[1])
  }

  def 'test mpx no query cross correlation false'() {
    given:
    var mp = new double[]{0.8106443110125229, 0.36123420440644344, 0.7520456848935667, 0.42833773532309266, 0.6871984848605489, 0.39204771374620745, 0.4028979513385712, 0.6541161237044746, 0.6155373168271888, 1.026609889502267, 0.8963078212971005, 0.6433281384283139, 0.5002237059826393, 0.8790165555463835, 0.4261574333144366, 0.7754136561837666, 0.9114091776576775, 0.4854349239569046, 0.4449744057351269, 1.0115056152957744, 1.659933435091076, 1.1947931467458364, 0.9333866174608633, 0.6433281384283139, 1.006827963280922, 1.0845268274947542, 0.8963078212971005, 0.4809499360572959, 0.5679781945824924, 1.0880723091374522, 1.265081184526383, 0.7329417256555751, 0.9600689711707685, 0.7425074741562939, 0.567531047582352, 0.515960623365906, 0.4493839590152485, 0.6625193162975406, 0.8478910670085921, 0.8906373784087319, 0.40957266604765424, 0.5353500319132934, 0.9759311981862917, 0.9670631085743358, 1.0880723091374522, 0.44627255649318603, 0.9183140029593616, 0.9558559329526902, 1.0321913309152035, 0.36123420440644344, 1.018944136007164, 0.40957266604765424, 0.5353500319132934, 0.5877069741794665, 0.4809499360572959, 0.42833773532309266, 0.6871984848605489, 0.39204771374620745, 0.5617691979076742, 0.6189084378964191, 0.9881871215933996, 0.7642636481956423, 0.4028979513385712, 0.6699514870145524, 0.582988609161628, 0.567531047582352, 0.6307890170937095, 0.6014757900601442, 0.9640127006187387, 1.0265294900568416, 0.4842839861959971, 0.8862122062229096, 0.894706312727349, 0.4854349239569046, 0.4449744057351269, 0.6014757900601442, 0.8906373784087319, 1.1203036986048183, 1.0824155404196334, 0.9550993928935435, 1.1115275279798638, 1.039351045869597, 1.5186274137076694, 1.1023607245964013, 0.515960623365906, 0.4493839590152485, 0.611473469822067, 0.7540834574554799, 0.992079864875988, 1.7261239601069096, 1.6098852660752319, 1.13652898187979, 0.992079864875988, 1.5673854937119218, 1.3208010566697452, 1.0321913309152035, 1.0265294900568416, 0.4842839861959971, 0.5002237059826393, 0.7157304302960894, 1.0115056152957744, 1.052020993230219, 0.53410772725371, 0.7157304302960894, 1.0170296626796977, 0.7329417256555751, 0.9600689711707685, 0.6763641538065726, 0.4261574333144366, 1.4078813621887736, 0.44627255649318603, 0.9640127006187387, 1.0090505359501154, 1.2769326813243422, 0.8936059628907048, 1.1356406334894906, 0.6085124794227836, 1.1811469738107099, 0.9670631085743358, 0.6189084378964191, 0.6085124794227836, 1.1018043642024244, 0.9759311981862917}
    var mpi = new int[] {52, 49, 54, 55, 56, 57, 62, 84, 85, 37, 26, 23, 98, 56, 108, 65, 72, 73, 74, 100, 101, 48, 53, 11, 47, 95, 10, 54, 55, 44, 104, 105, 106, 107, 65, 84, 85, 74, 75, 76, 51, 52, 122, 118, 29, 110, 23, 40, 95, 1, 46, 40, 41, 49, 27, 3, 4, 5, 6, 119, 64, 34, 6, 7, 14, 34, 74, 75, 111, 96, 97, 98, 33, 17, 18, 67, 39, 114, 102, 103, 10, 23, 55, 72, 35, 36, 65, 66, 92, 93, 71, 99, 88, 23, 21, 48, 69, 70, 12, 103, 19, 31, 98, 99, 100, 31, 32, 63, 14, 44, 45, 68, 1, 70, 3, 4, 120, 61, 43, 59, 116, 38, 42}

    when:
    var mpx = MPX.mpx(timeSeries, 6, false)

    then:
    equals(mpx.mp(), mp)
    equals(mpx.mpi(), mpi)
  }

  def 'test mpx no query cross correlation true'() {
    given:
    var mp = new double[]{0.945237983418586, 0.989125820797237, 0.9528689406527472, 0.984710565374857, 0.9606465202004472, 0.9871915491788643, 0.9864727700672652, 0.964344341392486, 0.9684261509660987, 0.9121726778980119, 0.9330526907901371, 0.9655107421921967, 0.9791480203310828, 0.9356108245896143, 0.9848658201692376, 0.949894471816977, 0.9307777759067797, 0.9803627445502462, 0.9834998148533892, 0.9147380325187597, 0.770385082588895, 0.8810391113740985, 0.9273991185287473, 0.9655107421921967, 0.9155247876962992, 0.901983463370347, 0.9330526907901371, 0.9807239299172069, 0.9731167308732344, 0.9013415541740244, 0.8666307997131103, 0.9552330355660856, 0.9231889642162585, 0.9540568875685034, 0.9731590425025065, 0.9778153862613221, 0.9831711714483151, 0.9634223462943866, 0.9400900615405859, 0.9338970883484351, 0.9860208526022181, 0.9761166952775363, 0.9206298580338891, 0.9220657453362119, 0.9013415541740244, 0.9834034004434197, 0.9297249493307295, 0.9238616196199285, 0.9112150880319584, 0.989125820797237, 0.9134794039747178, 0.9860208526022181, 0.9761166952775363, 0.971216709375068, 0.9807239299172069, 0.984710565374857, 0.9606465202004472, 0.9871915491788643, 0.9737012806901807, 0.9680793621250512, 0.918623851059746, 0.9513250896705573, 0.9864727700672652, 0.9625970837539158, 0.9716770234656492, 0.9731590425025065, 0.9668421013261627, 0.9698522394976271, 0.9225566260871472, 0.9121864338369701, 0.9804557517261763, 0.9345523271284603, 0.9332917178304859, 0.9803627445502462, 0.9834998148533892, 0.9698522394976271, 0.9338970883484351, 0.8954099685743637, 0.9023647164881727, 0.9239820958078654, 0.8970422128785811, 0.9099791169541479, 0.807814231527963, 0.8987334027389414, 0.9778153862613221, 0.9831711714483151, 0.9688416829753135, 0.9526131782659991, 0.9179814618089701, 0.7517080061954033, 0.7840224525061567, 0.8923584894456074, 0.9179814618089701, 0.7952752261751196, 0.8546237140583404, 0.9112150880319584, 0.9121864338369701, 0.9804557517261763, 0.9791480203310828, 0.9573108292623479, 0.9147380325187597, 0.9077709858169086, 0.976227411307323, 0.9573108292623479, 0.9138042221024684, 0.9552330355660856, 0.9231889642162585, 0.9618776276204599, 0.9848658201692376, 0.834822505833457, 0.9834034004434197, 0.9225566260871472, 0.9151514179915654, 0.8641202439471521, 0.9334556985905147, 0.8925266959639657, 0.969142713532228, 0.8837409855214835, 0.9220657453362119, 0.9680793621250512, 0.969142713532228, 0.8988355952520409, 0.9206298580338891}
    var mpi = new int[] {52, 49, 54, 55, 56, 57, 62, 84, 85, 37, 26, 23, 98, 56, 108, 65, 72, 73, 74, 100, 101, 48, 53, 11, 47, 95, 10, 54, 55, 44, 104, 105, 106, 107, 65, 84, 85, 74, 75, 76, 51, 52, 122, 118, 29, 110, 23, 40, 95, 1, 46, 40, 41, 49, 27, 3, 4, 5, 6, 119, 64, 34, 6, 7, 14, 34, 74, 75, 111, 96, 97, 98, 33, 17, 18, 67, 39, 114, 102, 103, 10, 23, 55, 72, 35, 36, 65, 66, 92, 93, 71, 99, 88, 23, 21, 48, 69, 70, 12, 103, 19, 31, 98, 99, 100, 31, 32, 63, 14, 44, 45, 68, 1, 70, 3, 4, 120, 61, 43, 59, 116, 38, 42}

    when:
    var mpx = MPX.mpx(timeSeries, 6, true)

    then:
    equals(mpx.mp(), mp)
    equals(mpx.mpi(), mpi)
  }

  def 'test mpx with query cross correlation false'() {
    given:
    var mp = new double[] {2.3668965327355505, 3.406013201994227, 2.903258498613403, 0.6971893604918976, 2.564482199238794, 4.1739525510420865, 3.6415334685242047, 3.593577627937353, 4.137081417670802, 4.457607139901716, 4.431555859764931, 3.2138410067775474, 0.6047023031090517, 3.104179903176288, 4.4704413299157135, 4.125780530323841, 3.7386771849938, 3.759127813545472, 4.043050494354302, 4.786387167956756, 3.8406642856176427, 2.6615627425970345, 3.52591075681746, 3.214588771573635, 1.1673105773258807, 2.9681630725325707, 4.231684110778988, 2.8078035727047905, 0.0, 3.1069506036555383, 4.513042682841902, 3.2187064353818027, 2.144885633847406, 3.9349560340272185, 4.2250836634586255, 3.821761785345031, 3.846206390701191, 4.4064420547591725, 4.790377123441244, 3.7430878242610315, 1.9099407632388654, 1.929431113798877, 3.4722046944834943, 2.3718192858194205, 3.1911690614398047, 4.053820330556801, 3.381534709024695, 1.9174055058898298, 3.1006973624748815, 3.6546095561285634, 2.9117004242206592, 1.6556301961839561, 1.9580427423249709, 3.6387515556921404, 3.1606124690559647, 0.5679781945824901, 2.643375755756869, 4.202408197533751, 3.7523411112426786, 3.551779992997021, 4.4967822474339245, 4.468366408151629, 3.4920148506393947, 3.6050715839334413, 4.341144702142613, 4.2001492696891916, 3.8578490765042703, 4.594175173745522, 3.712677713373187, 3.230294828773182, 1.9189154504498602, 1.6834692361706658, 3.4603195050456343, 4.021019724347909, 4.096874168432388, 4.605768013107246, 4.156449059287271, 2.0402285396416087, 1.1629816984710768, 4.04282959321277, 4.6316191983578765, 3.0044301207349906, 1.6435869891044308, 2.873140463889951, 3.8757666746585686, 4.028357105034513, 4.207940924791518, 4.002642512635051, 4.663606094977323, 3.719173514908543, 2.6233748670243293, 2.8037366924402285, 4.738290579046348, 3.0049182185796113, 2.435347999887961, 3.1701103260299917, 3.8572157511967204, 2.2161492289998006, 0.9865620496358983, 3.234225731274196, 4.873465031798678, 3.516413788704127, 1.175217809168927, 3.578334129818575, 4.7079896222663615, 3.271082546549932, 2.7586795054885647, 4.001399279436594, 4.500781547441157, 3.929933232627963, 3.842705957813382, 3.2338408014374944, 2.765570998878233, 2.36971926627377, 1.278356113689702, 1.6938320830183713, 4.191468434347383, 4.417641498556871, 3.1458027081346582, 3.396909990444241, 4.3130770219278975, 4.804077187095011, 3.5059890582543676}
    var mpi = new int[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

    when:
    var mpx = MPX.mpx(timeSeries, query, 6, false)

    then:
    equals(mpx.mp(), mp)
    equals(mpx.mpi(), mpi)
  }

  def 'test mpx with query cross correlation true'() {
    given:
    var mp = new double[] {0.5331500669437024, 0.03325617232008617, 0.2975908408524209, 0.9594939163014082, 0.45195258748227995, -0.45182332486256205, -0.105063833531827, -0.07615001400098792, -0.42628688803642134, -0.655855117808563, -0.6365572781847415, 0.13926883192959005, 0.9695279270512174, 0.19700559405970447, -0.6654038070182148, -0.4185054153666057, -0.1648089244660971, -0.17758682654759708, -0.3621881083248805, -0.9091251767984242, -0.2292251795682398, 0.4096736472682793, -0.03600388875342295, 0.1388682524727257, 0.8864488346719266, 0.26583399790450074, -0.4922625344516133, 0.3430199247588511, 1.0, 0.19557149553704048, -0.6972961880960691, 0.13666074023598082, 0.6166221348095843, -0.2903232491439346, -0.48761099693541315, -0.2171552619936366, -0.23277529998922364, -0.6180609651625197, -0.9123094153991006, -0.1675588716775987, 0.6960105234098783, 0.6897746314254021, -0.004683786699434654, 0.5312061062845879, 0.15137000177578314, -0.36945493936963725, 0.047101917638439024, 0.6936296771652805, 0.1988063221951094, -0.11301425064551801, 0.2935000532994361, 0.7715740544569896, 0.6805057182690423, -0.10337607367099752, 0.16754406837066327, 0.9731167308732346, 0.4177137178230668, -0.47168622155823936, -0.17333865126016146, -0.0512617598878265, -0.6850875484030747, -0.6638581964581578, -0.016180643090506076, -0.0830450937736976, -0.5704614437450731, -0.4701044906392209, -0.24024995809040417, -0.7588704605883081, -0.14866465028149628, 0.13043294326676988, 0.6931469578354008, 0.7638276109055796, 0.002182410250061261, -0.34738330196624406, -0.3986981626640471, -0.7677582492134889, -0.4396723985375035, 0.6531222921693224, 0.8872894640851108, -0.3620392599797441, -0.7876580332164381, 0.24778330413502747, 0.7748851507705526, 0.31208865622985293, -0.2517972763661616, -0.35230508047350334, -0.4755639022112746, -0.33509559032945274, -0.8124351507591361, -0.15268763616643038, 0.4264920255887569, 0.34492171328869403, -0.8709498009566146, 0.24753887497069463, 0.5057566766201422, 0.16253337673316826, -0.23984277927333986, 0.5907235495669658, 0.918891276851518, 0.12831532659699083, -0.979221784680374, -0.03043049444904258, 0.8849052584176823, -0.06703959538537155, -0.8470971902806467, 0.10833491447136795, 0.3658072821664473, -0.3342663494563077, -0.6880862114822347, -0.2870312677428058, -0.2305324231845384, 0.12852280591317525, 0.36263475418030444, 0.532035883254242, 0.8638171372160135, 0.7609110728781371, -0.4640339696775418, -0.6262963674809827, 0.1753271101243876, 0.03841687640167552, -0.5502194497568685, -0.9232631349638927, -0.024329939716612035}
    var mpi = new int[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

    when:
    var mpx = MPX.mpx(timeSeries, query, 6, true)

    then:
    equals(mpx.mp(), mp)
    equals(mpx.mpi(), mpi)
  }
}
