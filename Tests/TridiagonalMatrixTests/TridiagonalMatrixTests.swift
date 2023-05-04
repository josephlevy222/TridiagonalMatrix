import XCTest
@testable import TridiagonalMatrix
import Numerics
final class TridiagonalMatrixTests: XCTestCase {
    func testNormal() throws {
        try testExample(Complex<Double>(2.0,0.0))
    }
    func testSingular() throws {
        try testExample(Complex<Double>(1.0,0.0))
    }
    func testNearlySingular() throws {
        try testExample(Complex<Double>(1.732051,0.0))
    }
    
    func testExample(_ d: Complex<Double> ) throws {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct
        // results.
     
        let one = Complex<Double>(1.0,0.0)
        //let d = Complex<Double>(1.0,0.0) // nearly singular
        let lower = [-one,-one,-one,-one]
        let upper = lower
        let diagonal = [d,d,d,d,d]
        let tridiag = TridiagonalMatrix(diagonal: diagonal, upper: upper, lower: lower)
        let tridiagLU = TridiagonalLUMatrix(tridiag)
        let i0 = [one,0,0,0,0]
        let i1 = [0,one,0,0,0]
        let i2 = [0,0,one,0,0]
        let i3 = [0,0,0,one,0]
        let i4 = [0,0,0,0,one]
        let x0 = tridiagLU.solve(i0)
        let x1 = tridiagLU.solve(i1)
        let x2 = tridiagLU.solve(i2)
        let x3 = tridiagLU.solve(i3)
        let x4 = tridiagLU.solve(i4)
        //print(tridiag)
        let ii0 = tridiag * x0; //print(ii0)
        let ii1 = tridiag * x1; //print(ii1)
        let ii2 = tridiag * x2; //print(ii2)
        let ii3 = tridiag * x3; //print(ii3)
        let ii4 = tridiag * x4; //print(ii4)
        let tiny = 1e-15/tridiagLU.smallestPivot
        print("tiny=\(tiny), smallestPivot=\(tridiagLU.smallestPivot)")
        let e0 =  zip(i0, ii0).map { $0 - $1 } // error (difference) in identity matrix
        let e1 =  zip(i1, ii1).map { $0 - $1 }
        let e2 =  zip(i2, ii2).map { $0 - $1 }
        let e3 =  zip(i3, ii3).map { $0 - $1 }
        let e4 =  zip(i4, ii4).map { $0 - $1 }
        var result = e0.reduce(true) { $0 && $1.magnitude < tiny }
        result = e1.reduce(result) { $0 && $1.magnitude < tiny }
        result = e2.reduce(result) { $0 && $1.magnitude < tiny }
        result = e3.reduce(result) { $0 && $1.magnitude < tiny }
        result = e4.reduce(result) { $0 && $1.magnitude < tiny }
        var maxError = e0.reduce(0) { max($0,$1.magnitude) }
        maxError = e1.reduce(maxError) { max($0,$1.magnitude) }
        maxError = e2.reduce(maxError) { max($0,$1.magnitude) }
        maxError = e3.reduce(maxError) { max($0,$1.magnitude) }
        maxError = e4.reduce(maxError) { max($0,$1.magnitude) }
        print("maxError=\(maxError)")
        let det = tridiagLU.determinate(); print("det=\(det)")
        switch d {
        case Complex<Double>(2.0,0.0) :
            XCTAssertEqual(Complex<Double>(6.0, 0.0), det)
        case  Complex<Double>(1.0,0.0) :
            XCTAssertEqual(Complex<Double>(0.0, 0.0), det)
        
        default :
            XCTAssertLessThan(.zero, det.magnitude)
        }
            
       XCTAssertTrue(result || tridiagLU.singular)
    }
}
