Function EigenValueJ(x As Object)
'Real Eigen Value by Jacobi method
'This function for Real symmetric matrix only
'original 6/14/2004 , modified 12/28/2004
'copy 2/18/2014

    '
    Dim i, j, k, N, II, n_div, K_new, i_div, i_rot, j_rot, Kmax As Integer
    Dim a(), U(), v() As Double
    Dim b1, b2, b3, a_max, EpsRate, Eps, K_iter, ax, v_cos, v_sin, v_tan As Double
    Dim i1, i2 As Double
    
    N = x.Columns.Count
    ReDim a(1 To N, 1 To N), v(1 To N, 1 To N), a2(1 To N)
    
    'define acceptable error rate
        Eps = 10 ^ -16
    'max itaretion number
        Kmax = 4000
        
    'initialize matrix A()
        For i = 1 To N
            For j = 1 To N
                a(i, j) = x.Cells(i, j).Value
            Next j
        Next i
        
    'initialize matrix V()
    'V == unit matrix
        For i = 1 To N
            For j = 1 To N
                v(i, j) = 0
            Next j
        Next i
        
        For i = 1 To N
            v(i, i) = 1
        Next i
        
    Do
        kiter = kiter + 1 'add itaration number
        knew = 0 'flag :if nondaigonal elemnt has unacceptable error, flag was 1
    
            ax = 0
            
            'looking for maximum element and replace ax
            For i = 1 To N
                For j = i + 1 To N
                    b1 = Abs(a(i, j))
                    If b1 > ax Then
                        irot = i
                        jrot = j
                        ax = b1
                    End If
                Next j
            Next i
            
            If ax > Eps Then
            knew = 1
            
            'calcurate cosine and sine
            b1 = a(irot, irot) - a(jrot, jrot)
            b2 = -b1 - Sqr(b1 ^ 2 + 4 * ax ^ 2)
            b3 = 2 * a(irot, jrot)
            vtan = b2 / b3
            vcos = 1 / Sqr(1 + vtan ^ 2)
            vsin = vcos * vtan
            
            'rotate V vector
            For i = 1 To N
                b1 = v(i, irot)
                v(i, irot) = vcos * b1 + vsin * v(i, jrot)
                v(i, jrot) = vcos * v(i, jrot) - vsin * b1
            Next i
            
            'rotate A matrix
            For i = 1 To irot - 1
                b1 = a(i, irot)
                a(i, irot) = vcos * b1 + vsin * a(i, jrot)
                a(i, jrot) = vcos * a(i, jrot) - vsin * b1
            Next i
            
            For i = irot + 1 To jrot - 1
                b1 = a(irot, i)
                a(irot, i) = vcos * b1 + vsin * a(i, jrot)
                a(i, jrot) = vcos * a(i, jrot) - vsin * b1
            Next i
            
            For i = jrot + 1 To N
                b1 = a(irot, i)
                a(irot, i) = vcos * b1 + vsin * a(jrot, i)
                a(jrot, i) = vcos * a(jrot, i) - vsin * b1
            Next i
            
            'erase nondiagonal element
            b1 = a(irot, irot)
            b2 = 2 * vcos * vsin * a(irot, jrot)
            b3 = vcos ^ 2 * b1 + b2
            
            a(irot, irot) = b3 + vsin ^ 2 * a(jrot, jrot)
            b3 = vcos ^ 2 * a(jrot, jrot) + vsin ^ 2 * b1
            a(jrot, jrot) = b3 - b2
            a(irot, jrot) = 0
            a(jrot, irot) = a(irot, jrot)
        End If
    Loop While kiter < Kmax And knew > 0
    
    'output
    For i = 1 To N
        a2(i) = a(i, i)
    Next
    
'Eigen values/vectors
EigenValueJ = a2 'if you need Eigen vectors ,replace a2 to v

End Function
