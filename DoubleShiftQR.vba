Function shiftQR(x As Object)
'Extraction Ak(Eigen values) by Double shift QR method after Hessenberg
'Note: This program is NOT stable and NOT smart. Need improvement for convergence.
'refer from `Numerical recipes in C`
'https://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf
'　6/10/2014
'wr real
'wi imaginary

Dim i, j, k, l, N, nn, M, its, mmin As Integer
Dim a(), wr(), wi(), eig() As Double
Dim U, r, q, p, y, w, z, xx, t, anorm As Double

N = x.Columns.Count
ReDim a(1 To N, 1 To N)
ReDim wr(1 To N), wi(1 To N), eig(1 To 2, 1 To N)

For i = 1 To N
    For j = 1 To N
        a(i, j) = x(i, j)
    Next j
Next i

'MsgBox (a(4, 3))

'Matrix norm include sub-diagonal element
anorm = 0#
For i = 1 To N
    If i - 1 = 0 Then
        For j = 1 To N
            anorm = anorm + Abs(a(i, j))
        Next j
    Else
        For j = i - 1 To N
            anorm = anorm + Abs(a(i, j))
        Next j
    End If
Next i

nn = N
t = 0#
'search for Eigen values
Do While nn >= 1
'MsgBox (nn)
    its = 0
        Do  'Iteration start
            'Search for extreme small value (almost 0) on sub-diagnal term
            For l = nn To 2 Step -1
                
                s = Abs(a(l - 1, l - 1)) + Abs(a(l, l))
                 If s = 0 Then
                    s = anorm
                End If
                If Abs(a(l, l - 1)) + s = s Then

                    a(l, l - 1) = 0#

                    Exit For
                End If
                
                If l = 2 Then
                    Exit For
                End If
            Next l
            
            xx = a(nn, nn)
            'Real Eigen value --- when a(n,n-1) = 0
            If l = nn Then
                wr(nn) = xx + t
                wi(nn) = 0
                nn = nn - 1
            Else
                y = a(nn - 1, nn - 1)
                w = a(nn, nn - 1) * a(nn - 1, nn)
                
                'Complex Eigen value --- extraction from 2 * 2 matrix (right of 0 sub-diagnal term)
                If l = (nn - 1) Then
                    p = 0.5 * (y - xx)
                    q = p * p + w
                    z = Sqr(Abs(q))
                    xx = xx + t
                    
                    'Real conjugation
                    If q >= 0# Then
                        z = p + Abs(z) * (p / Abs(p))
                        wr(nn - 1) = xx + z
                        wr(nn) = wr(nn - 1)
                        If z <> 0# Then
                            wr(nn) = xx - w / z
                        End If
                        wi(nn - 1) = 0#
                        wi(nn) = 0#
                    Else 'Complex conjugation
                        wr(nn - 1) = xx + p
                        wr(nn) = wr(nn - 1)
                        wi(nn - 1) = -z
                        wi(nn) = z
                    End If
                    
                    nn = nn - 2
                    
                    Exit Do
                Else 'Continue iteration when Eigen value is not extracted
                    If its = 30 Then
                        MsgBox ("Too many iterations")
                    End If
                    If its = 10 Or its = 20 Then 'form exceptional shift
                        t = t + xx
                        For i = 1 To nn
                            a(i, i) = a(i, i) - xx
                        Next i
                        s = Abs(a(nn, nn - 1)) + Abs(a(nn - 1, nn - 2))
                        y = 0.75 * s
                        xx = y
                        w = -0.4375 * s * s
                    End If
                    
                    its = its + 1
                    
                    'Search extreme small value with form exeptional shift and sub-diagonal term
                    For M = nn - 2 To l Step -1
                        z = a(M, M)
                        r = xx - z
                        s = y - z
                        p = (r * s - w) / a(M + 1, M) + a(M, M + 1) 'formula 11.6.23
                        q = a(M + 1, M + 1) - z - r - s
                        r = a(M + 2, M + 1)
                        s = Abs(p) + Abs(q) + Abs(r) 'Scaling to prevent to overflow and underflow
                        p = p / s
                        q = q / s
                        r = r / s
                        If M = l Then
                            Exit For
                        End If
                        
                        U = Abs(a(M, M - 1)) * (Abs(q) + Abs(r))
                        v = Abs(p) * (Abs(a(M - 1, M - 1)) + Abs(z) + Abs(a(M + 1, M + 1)))
                        'formula 11.6.26
                        
                        If U + v = v Then
                            Exit For
                        End If
                        
                            
                    Next M
                    
                    For i = M + 2 To nn
                        a(i, i - 2) = 0#
                        If i <> M + 2 Then
                            a(i, i - 3) = 0#
                        End If
                    Next i
                    
                    For k = M To nn - 1  'do?
                        'Double QR step row L~nn and row m~nn
                        If k <> M Then
                            'Vector for Householder
                            p = a(k, k - 1)
                            q = a(k + 1, k - 1)
                            r = 0#
                            If k <> (nn - 1) Then
                                r = a(k + 2, k - 1)
                            End If
                            
                                xx = Abs(p) + Abs(q) + Abs(r)
                            If xx <> 0# Then
                                'Scaling to prevent to overflow and underflow

                                p = p / xx
                                q = q / xx
                                r = r / xx
                                
                            End If
                        End If
                            s = Sqr(p * p + q * q + r * r) * (p / Abs(p))
                        If s <> 0# Then

                            If k = M Then
                                If l <> M Then
                                    a(k, k - 1) = -a(k, k - 1)
                                Else
                                    a(k, k - 1) = -s * xx
                                
                                p = p + s 'formula 11.6.24
                                xx = p / s
                                y = q / s
                                z = r / s
                                q = q / p
                                r = r / p
                                
                                'transform row
                                For j = k To nn
                                    p = a(k, j) + q * a(k + 1, j)
                                    If k <> (nn - 1) Then
                                        p = p + r * a(k + 2, j)
                                        a(k + 2, j) = a(k + 2, j) - p * z
                                    End If
                                    
                                    a(k + 1, j) = a(k + 1, j) - p * y
                                    a(k, j) = a(k, j) - p * xx
                                Next j

                                    If nn < k + 3 Then
                                        mmin = nn
                                    Else
                                        mmin = k + 3
                                    End If
                                
                                'transform column
                                For i = l To mmin
                                    p = xx * a(i, k) + y * a(i, k + 1)
                                    If k <> (nn - 1) Then
                                        p = p + z * a(i, k + 2)
                                        a(i, k + 2) = a(i, k + 2) - p * r
                                    End If
                                    
                                    a(i, k + 1) = a(i, k + 1) - p * q
                                    a(i, k) = a(i, k) - p
                                Next i
                                
                                End If
                            End If
                        End If
                    Next k
                End If
            End If
        Loop While l < nn - 1
Loop
                                
For i = 1 To N
    eig(1, i) = wr(i)
    eig(2, i) = wi(i)
Next i

    shiftQR = a
End Function
