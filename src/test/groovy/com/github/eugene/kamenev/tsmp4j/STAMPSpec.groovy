package com.github.eugene.kamenev.tsmp4j

class STAMPSpec extends BaseSpec {

  def 'test stamp with no query'() {
    given:
    var mp = new double[] {0.8106443, 0.3612342, 0.7520457, 0.4283377, 0.6871985, 0.3920477, 0.4028980, 0.6541161, 0.6155373, 1.0266099, 0.8963078, 0.6433281, 0.5002237, 0.8790166, 0.4261574, 0.7754137, 0.9114092, 0.4854349, 0.4449744, 1.0115056, 1.6599334, 1.1947931, 0.9333866, 0.6433281, 1.0068280, 1.0845268, 0.8963078, 0.4809499, 0.5679782, 1.0880723, 1.2650812, 0.7329417, 0.9600690, 0.7425075, 0.5675310, 0.5159606, 0.4493840, 0.6625193, 0.8478911, 0.8906374, 0.4095727, 0.5353500, 0.9759312, 0.9670631, 1.0880723, 0.4462726, 0.9183140, 0.9558559, 1.0321913, 0.3612342, 1.0189441, 0.4095727, 0.5353500, 0.5877070, 0.4809499, 0.4283377, 0.6871985, 0.3920477, 0.5617692, 0.6189084, 0.9881871, 0.7642636, 0.4028980, 0.6699515, 0.5829886, 0.5675310, 0.6307890, 0.6014758, 0.9640127, 1.0265295, 0.4842840, 0.8862122, 0.8947063, 0.4854349, 0.4449744, 0.6014758, 0.8906374, 1.1203037, 1.0824155, 0.9550994, 1.1115275, 1.0393510, 1.5186274, 1.1023607, 0.5159606, 0.4493840, 0.6114735, 0.7540835, 0.9920799, 1.7261240, 1.6098853, 1.1365290, 0.9920799, 1.5673855, 1.3208011, 1.0321913, 1.0265295, 0.4842840, 0.5002237, 0.7157304, 1.0115056, 1.0520210, 0.5341077, 0.7157304, 1.0170297, 0.7329417, 0.9600690, 0.6763642, 0.4261574, 1.4078814, 0.4462726, 0.9640127, 1.0090505, 1.2769327, 0.8936060, 1.1356406, 0.6085125, 1.1811470, 0.9670631, 0.6189084, 0.6085125, 1.1018044, 0.9759312}
    var mpi = new int[] {52, 49, 54, 55, 56, 57, 62, 84, 85, 37, 26, 23, 98, 56, 108, 65, 72, 73, 74, 100, 101, 48, 53, 11, 47, 95, 10, 54, 55, 44, 104, 105, 106, 107, 65, 84, 85, 74, 75, 76, 51, 52, 122, 118, 29, 110, 23, 40, 95, 1, 46, 40, 41, 49, 27, 3, 4, 5, 6, 119, 64, 34, 6, 7, 14, 34, 74, 75, 111, 96, 97, 98, 33, 17, 18, 67, 39, 114, 102, 103, 10, 23, 55, 72, 35, 36, 65, 66, 92, 93, 71, 99, 88, 23, 21, 48, 69, 70, 12, 103, 19, 31, 98, 99, 100, 31, 32, 63, 14, 44, 45, 68, 1, 70, 3, 4, 120, 61, 43, 59, 116, 38, 42}

    when:
    var stamp = STAMP.stamp(timeSeries, 6)

    then:
    equals(stamp.distances(), mp)
    equals(stamp.indexes(), mpi)
  }

  def 'test stamp with query'() {
    given:
    var mp = new double[] {2.3668965327355567, 3.4060132019942193, 2.9032584986134187, 0.697189360491946, 2.564482199238838, 4.173952551042108, 3.641533468524242, 3.5935776279373743, 4.137081417670844, 4.457607139901744, 4.431555859764974, 3.213841006777569, 0.6047023031090474, 3.1041799031762967, 4.470441329915717, 4.125780530323814, 3.738677184993776, 3.759127813545457, 4.043050494354343, 4.786387167956826, 3.8406642856176942, 2.661562742597054, 3.525910756817523, 3.214588771573668, 1.1673105773258903, 2.9681630725325654, 4.231684110779005, 2.8078035727048456, 6.542950391799256E-7, 3.1069506036555286, 4.513042682841856, 3.218706435381751, 2.14488563384738, 3.934956034027165, 4.225083663458599, 3.821761785345017, 3.8462063907011657, 4.406442054759158, 4.790377123441251, 3.743087824261052, 1.909940763238877, 1.9294311137988711, 3.4722046944838154, 2.3718192858198504, 3.1911690614398625, 4.05382033055673, 3.3815347090246832, 1.9174055058897412, 3.1006973624748007, 3.6546095561285057, 2.9117004242206312, 1.6556301961839521, 1.9580427423249531, 3.6387515556921874, 3.1606124690560313, 0.5679781945823819, 2.6433757557568747, 4.202408197533685, 3.7523411112426395, 3.551779992997057, 4.496782247434107, 4.468366408151725, 3.4920148506394533, 3.6050715839335767, 4.3411447021426675, 4.200149269689229, 3.8578490765041473, 4.594175173745556, 3.712677713373182, 3.2302948287732107, 1.9189154504499335, 1.683469236170608, 3.460319505045696, 4.021019724347921, 4.0968741684324295, 4.605768013107236, 4.156449059287281, 2.04022853964161, 1.162981698471108, 4.042829593212975, 4.631619198358055, 3.0044301207350115, 1.6435869891044623, 2.8731404638900333, 3.8757666746586934, 4.028357105034631, 4.2079409247915605, 4.002642512635227, 4.6636060949775064, 3.7191735149087055, 2.6233748670246717, 2.8037366924403546, 4.738290579046625, 3.0049182185798626, 2.435347999888012, 3.1701103260300556, 3.857215751196695, 2.216149228999814, 0.9865620496359858, 3.234225731274193, 4.87346503179918, 3.5164137887045586, 1.1752178091695555, 3.5783341298193423, 4.707989622267232, 3.2710825465503457, 2.758679505488824, 4.001399279436573, 4.5007815474411155, 3.929933232627954, 3.842705957813346, 3.233840801437425, 2.7655709988782506, 2.369719266273785, 1.2783561136897605, 1.6938320830185023, 4.191468434347509, 4.417641498556896, 3.1458027081346565, 3.396909990444282, 4.313077021928001, 4.804077187095176, 3.5059890582542703}
    var mpi = new int[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

    when:
    var stamp = STAMP.stamp(timeSeries, query, 6)

    then:
    equals(stamp.distances(), mp)
    equals(stamp.indexes(), mpi)
  }
}
