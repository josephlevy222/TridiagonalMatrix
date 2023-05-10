//
//  TridiagonalMatrix.swift
//
//  Created by Joseph Levy on 4/19/23.
//

import Foundation

import Numerics
// Set up protocol for Real and Complex<Real> types
protocol RealOrComplex : AlgebraicField where Magnitude : Real {}
// and enumerate types that fit it that you want to use
extension Complex : RealOrComplex {}
extension Float   : RealOrComplex {}
extension Double  : RealOrComplex {}

struct TridiagonalMatrix<T: RealOrComplex > {
    var lower: [T]
    var diagonal: [T]
    var upper: [T]
    let size: Int
    init(diagonal: [T], upper: [T], lower: [T]) {
        precondition(diagonal.count > 0, "Diagonal must not be empty")
        precondition(diagonal.count == upper.count + 1, "Invalid upper size")
        precondition(diagonal.count == lower.count + 1, "Invalid lower size")
        self.diagonal = diagonal
        self.upper = upper
        self.lower = lower
        size = diagonal.count
    }
}

typealias ColumnVector<T: RealOrComplex> = Array<T>

func *<T: RealOrComplex>(_ A: TridiagonalMatrix<T>, _ x: ColumnVector<T>) -> ColumnVector<T> {
    precondition(x.count == A.size, "Invalid column vector size")
    let n = x.count
    var b : ColumnVector<T> = x
    if n==0 { return [] }
    b[0] = A.diagonal[0]*x[0]+A.upper[0]*x[1]
    if n==1 { return b }
    b[n-1] = A.lower[n-2]*x[n-2]+A.diagonal[n-1]*x[n-1]
    if n==2 { return b }
    for j in 1..<n-1 {
        b[j] = A.lower[j-1]*x[j-1] + A.diagonal[j]*x[j] + A.upper[j]*x[j+1]
    }
    return b
}
/// implements Ax+y where A is a TridiagonalMatrix and x and y are ColumnVector types
func AXpY<T: RealOrComplex>(A: TridiagonalMatrix<T>, x: ColumnVector<T>, y: ColumnVector<T>) -> ColumnVector<T> {
    precondition(x.count == A.size, "Invalid x vector size")
    precondition(y.count == A.size, "Invalid y vector size")
    let n = x.count
    if n==0 { return [] }
    var b = y
    b[0] += A.diagonal[0]*x[0]+A.upper[0]*x[1]
    if n==1 { return b }
    b[n-1] += A.lower[n-2]*x[n-2]+A.diagonal[n-1]*x[n-1]
    if n==2 { return b }
    for j in 1..<n-1 {
        b[j] += A.lower[j-1]*x[j-1] + A.diagonal[j]*x[j] + A.upper[j]*x[j+1]
    }
    return b
}

func aXpY<T: RealOrComplex >(a: T , x: ColumnVector<T>, y: ColumnVector<T>) -> ColumnVector<T> {
    precondition(x.count == y.count, "Vector size mismatchector")
    return y + x.map { a*$0 }
}

struct TridiagonalLUMatrix<T: RealOrComplex > {
    private var au0 = [T]()
    private var au1 = [T]()
    private var au2 = [T]()
    private var al  = [T]()
    private var indx = [Int]()
    private var d : Bool = true
    var approximateConditionNumber : T.Magnitude = .infinity
    init(_ A: TridiagonalMatrix<T>) {
        let n = A.size
        if n==0 { return }
        al = Array(repeating: T.zero, count: n)
        indx = Array(0..<n)
        au0 = [A.diagonal[0]] + A.lower
        au1 = A.diagonal
        au2 = A.upper + [0]
        au1[0] = au2[0]
        au2[0] = T.zero // completes rearrangement
        let maxElement = [au0,au1,au2].map {$0.map {$0.magnitude}.max()!}.max()!
        debugPrint("maxElement=\(maxElement)")
        approximateConditionNumber = au0[0].reciprocal?.magnitude ?? .infinity
        for k in 0..<n-1 {
            // Find the biggest pivot and check if 0
            var au0max = au0[k].magnitude
            var i = k
            let au0p1mag = au0[k+1].magnitude
            if au0p1mag > au0max {
                au0max = au0p1mag
                i = k+1
            }
            if let a =  au0[i].reciprocal?.magnitude,
               approximateConditionNumber < a { approximateConditionNumber = a }
            indx[k] = i
            if i != k {
                d.toggle() // to get sign of determinate
                au0.swapAt(k, i)
                au1.swapAt(k, i)
                au2.swapAt(k, i)
            }
            
            al[k] = au0[k+1]/au0[k] // fix for singular (maybe?)
            au0[k+1] = au1[k+1] - al[k] * au1[k]
            au1[k+1] = au2[k+1] - al[k] * au2[k]
            au2[k+1] = T.zero
            
        }
        indx[n-1] = n-1
        if let a=au0[n-1].reciprocal?.magnitude, approximateConditionNumber < a { approximateConditionNumber = a}
        approximateConditionNumber *= maxElement
    }
    
    func solveInPlace(_ x: inout ColumnVector<T>) {
        let n = au0.count
        if n==0 {return}
        for k in 0..<n {
            x.swapAt(k, indx[k])
            if k<n-1 { x[k+1] -= al[k]*x[k] }
        }
        x[n-1] = x[n-1]/au0[n-1]
        if n==1 {return}
        x[n-2] = (x[n-2]-au1[n-2]*x[n-1])/au0[n-2]
        if n==2 {return}
        for i in (0..<n-2).reversed() {
            x[i] = (x[i]-au1[i]*x[i+1]-au2[i]*x[i+2])/au0[i]
        }
        return
    }
    
    func solve(_ b: ColumnVector<T>) -> ColumnVector<T> {
        var x = b
        solveInPlace(&x)
        return x
    }
    
    func determinate() -> T {
        var det = T(exactly: 1)!
        if !d { det = -det}
        for i in 0..<au0.count {
            det *= au0[i]
        }
        return det
    }
}
