//
//  TridiagonalMatrix.swift
//
//  Created by Joseph Levy on 4/19/23.
//

import Foundation

import Numerics
struct TridiagonalMatrix<T: AlgebraicField > {
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

typealias ColumnVector<T> = Array<T>

func *<T: AlgebraicField>(_ A: TridiagonalMatrix<T>, _ x: ColumnVector<T>) -> ColumnVector<T> {
    precondition(x.count == A.size, "Invalid column vector size")
    var b : ColumnVector<T> = x
    let n = x.count
    b[0] = A.diagonal[0]*x[0]+A.upper[0]*x[1]
    for j in 1..<n-1 {
        b[j] = A.lower[j-1]*x[j-1] + A.diagonal[j]*x[j] + A.upper[j]*x[j+1]
    }
    b[n-1] = A.lower[n-2]*x[n-2]+A.diagonal[n-1]*x[n-1]
    return b
}
/// implements Ax+y where A is a TridiagonalMatrix and x and y are ColumnVector types
func aXpY<T: AlgebraicField>(A: TridiagonalMatrix<T>, x: ColumnVector<T>, y: ColumnVector<T>) -> ColumnVector<T> {
    precondition(x.count == A.size, "Invalid x vector size")
    precondition(y.count == A.size, "Invalid y vector size")
    var b : ColumnVector<T> = y
    let n = x.count
    b[0] += A.diagonal[0]*x[0]+A.upper[0]*x[1]
    for j in 1..<n-1 {
        b[j] += A.lower[j-1]*x[j-1] + A.diagonal[j]*x[j] + A.upper[j]*x[j+1]
    }
    b[n-1] += A.lower[n-2]*x[n-2]+A.diagonal[n-1]*x[n-1]
    return b
}

struct TridiagonalLUMatrix<T: AlgebraicField> {
    private var au0 : [T]
    private var au1 : [T]
    private var au2 : [T]
    private var al  : [T]
    private var indx : [Int]
    private var d : Bool
    var singular : Bool
    var smallestPivot : T.Magnitude
    init(_ A: TridiagonalMatrix<T>) {
        let n = A.size
        al = Array(repeating: T.zero, count: n)
        indx = Array(0..<n)
        d = true
        singular = false
        au0 = [A.diagonal[0]] + A.lower
        au1 = A.diagonal
        au2 = A.upper + [0]
        au1[0] = au2[0]
        au2[0] = T.zero // completes rearrangement
        smallestPivot = au0[0].magnitude
        for k in 0..<n {
            // Find the biggest pivot and check if 0
            var au0max = au0[k].magnitude
            var i = k
            if k<n-1 { //for j in k+1..<el { // j = k+1 or none
                let au0p1mag = au0[k+1].magnitude
                if au0p1mag > au0max {
                    au0max = au0p1mag
                    i = k+1
                }
            }
            if smallestPivot > au0max { smallestPivot = au0max }
            indx[k] = i+1
            if au0max == 0 { singular = true; print("singular") }  // want au0[k]=tinyValue
            if i != k {
                d.toggle() // to get sign of determinate
                au0.swapAt(k, i)
                au1.swapAt(k, i)
                au2.swapAt(k, i)
            }
            if k<n-1 { //for i in k+1..<el { // i is k+1 el is k+1 or k if k+1 is n-1
                al[k] = au0[k+1]/au0[k] // fix for singular (maybe?)
                au0[k+1] = au1[k+1] - al[k] * au1[k]
                au1[k+1] = au2[k+1] - al[k] * au2[k]
                au2[k+1] = T.zero
            }
        }
    }
    
    func solveInPlace(_ x: inout ColumnVector<T>) {
        let n = au0.count
        for k in 0..<n {
            x.swapAt(k, indx[k]-1)
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
