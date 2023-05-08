import XCTest
@testable import TridiagonalMatrix
import Numerics

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
    
    func testExample<T : RealOrComplex >(_ d: T, det: T) throws where T.Magnitude : Real {
        guard let one = T(exactly: 1) else { print("Test of type like \(d) not possible"); return }
        let lower = [-one,-one,-one,-one]
        let upper = lower
        let diagonal = [d,d,d,d,d]
        let tridiag = TridiagonalMatrix(diagonal: diagonal, upper: upper, lower: lower)
        let tridiagLU = TridiagonalLUMatrix(tridiag)
        let i = [ [one,0,0,0,0],
                  [0,one,0,0,0],
                  [0,0,one,0,0],
                  [0,0,0,one,0],
                  [0,0,0,0,one] ]
        let x = i.map {icol in tridiagLU.solve(icol)}
        let ii = x.map {xcol in tridiag*xcol}
        
        var tolerance = T.Magnitude.ulpOfOne*2
        if let condition = tridiagLU.approximateConditionNumber {
            print("conditionNumber=\(condition)")
            tolerance *= condition
        }
        print("tolerance=\(tolerance)")
        let e = zip(i,ii).map { zip($0,$1).map { $0 - $1 } }

        let okay = e.flatMap { $0 }.reduce(true) { $0 && $1.magnitude < tolerance }

        let maxError = e.flatMap {$0}.reduce(T.Magnitude.zero) { max($0,$1.magnitude)}
        print("maxError=\(maxError)")
        
        let determinate = tridiagLU.determinate()
        print("determinate=\(determinate) vs. \(det)")
    
        XCTAssertTrue(determinate.isApproximatelyEqual(to: det))
        
        XCTAssertTrue(okay || tridiagLU.singular)
    }
}
