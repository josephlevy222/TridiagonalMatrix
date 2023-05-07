import XCTest
@testable import TridiagonalMatrix
import Numerics
func roundoff<T: RealOrComplex>(_ d: T) -> T.Magnitude {
    if d.magnitude is Double { return Double.ulpOfOne*2 as? T.Magnitude ?? .zero }
    if d.magnitude is Float { return Float.ulpOfOne*2 as? T.Magnitude ?? .zero }
    if d.magnitude is Float80 { return Float80.ulpOfOne*2 as? T.Magnitude ?? .zero}
    return Double.ulpOfOne*2 as? T.Magnitude ?? .zero
}

func aboutEquals<T: RealOrComplex>(_ d1: T, _ d2: T) -> Bool {
    if d1 is Complex<Double> { return (d1 as! Complex<Double>).isApproximatelyEqual(to: d2 as! Complex<Double>) }
    if d1 is Complex<Float> { return (d1 as! Complex<Float>).isApproximatelyEqual(to: d2 as! Complex<Float>)}
    if d1 is Complex<Float80> { return(d1 as! Complex<Float80>).isApproximatelyEqual(to: d2 as! Complex<Float80>)}
    if d1 is Float { return (d1 as! Float).isApproximatelyEqual(to: d2 as! Float)}
    if d1 is Double { return (d1 as! Double).isApproximatelyEqual(to: d2 as! Double)}
    if d1 is Float80 { return (d1 as! Float80).isApproximatelyEqual(to: d2 as! Float80)}
    return true
}

final class TridiagonalMatrixTests: XCTestCase {
    func testNormal() throws {
        try testExample(Complex<Double>(2.0,0.0), det: Complex(6.0))
    }
    func testSingular() throws {
        try testExample(Complex<Float>(1.0,0.0), det: 0)
    }
    func testNearlySingular() throws {
        try testExample(Complex<Float>(1.732051,0.0), det: Complex(2.265e-6))
    }
    
    func testExample<T : RealOrComplex >(_ d: T, det: T) throws {
        guard let one = T(exactly: 1) else { print("Test of type like \(d) not possible"); return }
        let lower = [-one,-one,-one,-one]
        let upper = lower
        let diagonal = [d,d,d,d,d]
        let tridiag = TridiagonalMatrix(diagonal: diagonal, upper: upper, lower: lower)
        let tridiagLU = TridiagonalLUMatrix(tridiag)
        var i = Array<[T]>(repeating: [], count: 5)
        i[0] = [one,0,0,0,0]
        i[1] = [0,one,0,0,0]
        i[2] = [0,0,one,0,0]
        i[3] = [0,0,0,one,0]
        i[4] = [0,0,0,0,one]
        var x = i
        for j in 0..<i.count {
            x[j] = tridiagLU.solve(i[j])
        }
        var ii = i
        for j in 0..<i.count {
            ii[j] = tridiag*x[j]
        }

        var tolerance : T.Magnitude = roundoff(d)
        if let condition = tridiagLU.approximateConditionNumber {
            print("conditionNumber=\(condition)")
            tolerance *= condition
        }
        print("tolerance=\(tolerance)")
        var e = i
        for j in 0..<i.count { // error (difference) in identity matrix
            e[j] = zip(i[j],ii[j]).map { $0 - $1 }
        }
        var result = true
        for j in 0..<i.count {
            result = e[j].reduce(result) { $0 && ($1.magnitude) < tolerance }
        }

        var maxError = T.Magnitude.zero
        for j in 0..<i.count {
            maxError = e[j].reduce(maxError) { max($0,$1.magnitude) }
        }
        print("maxError=\(maxError)")
        let determinate = tridiagLU.determinate()
        print("determinate=\(determinate) vs. \(det)")
    
        XCTAssertTrue(aboutEquals(determinate,det))
        
        XCTAssertTrue(result || tridiagLU.singular)
    }
}
