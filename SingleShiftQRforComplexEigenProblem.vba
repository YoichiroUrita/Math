Function singleQR(x As Object)
'Extraction Ak by Single shift QR method after Hessenberg
'pukiwiki for PBCG Lab 
'http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?%B8%C7%CD%AD%C3%CD/%B8%C7%CD%AD%A5%D9%A5%AF%A5%C8%A5%EB
'Y.Urita
'　6/10/2014

Dim i, j, k, N, cnt As Integer
Dim a(), r(), q(), t(), U(), v() As Double
Dim c, s, b1, rq, e, Eps As Double

N = x.Columns.Count
ReDim a(1 To N, 1 To N), r(1 To N, 1 To N), q(1 To N, 1 To N), t(1 To N, 1 To N)
ReDim U(1 To N), v(1 To N)

'Substitute Hessenberg matrix to R and A
For i = 1 To N
    For j = 1 To N
        r(i, j) = x(i, j)
        a(i, j) = x(i, j)
    Next j
Next i

Eps = 0.00000001
cnt = 0
itermax = 500 'max iteration number

Do While cnt < itermax
    'unit matrixization of Q
    For i = 1 To N
        For j = 1 To N
            If i = j Then
                q(i, j) = 1#
            Else
                q(i, j) = 0#
            End If
        Next j
    Next i
    
    For k = 1 To N - 1
        'sin,cos
        b1 = Sqr(r(k, k) * r(k, k) + r(k + 1, k) * r(k + 1, k))
        
        c = r(k, k) / b1
        s = -r(k + 1, k) / b1
        
    
         'calculate R matrix
        For j = k + 1 To N
            U(j) = c * r(k, j) - s * r(k + 1, j)
            v(j) = s * r(k, j) + c * r(k + 1, j)
        Next j
    
        r(k, k) = b1
        r(k + 1, k) = 0#
        For j = k + 1 To N
            r(k, j) = U(j)
            r(k + 1, j) = v(j)
        Next j
    
        'calculate Q matrix
        For j = 1 To k
            U(j) = c * q(k, j)
            v(j) = s * q(k, j)
        Next j
        q(k, k + 1) = -s
        q(k + 1, k + 1) = c
        For j = 1 To k
            q(k, j) = U(j)
            q(k + 1, j) = v(j)
        Next j
    Next k
    
    'A=RQ
    For i = 1 To N
        For j = 1 To N
            rq = 0#
            For M = 1 To N
                rq = rq + r(i, M) * q(j, M)
            Next M
            
            t(i, j) = rq
        Next j
    Next i
    
    'convergence check
    e = 0#
    For i = 1 To N
        e = e + Abs(t(i, i) - a(i, i))
    Next i
    
    If e < Eps Then
        'Exit Do
    End If
    
    'Substitute RQ to R,A
    For i = 1 To N
        For j = 1 To N
            r(i, j) = t(i, j)
            a(i, j) = t(i, j)
        Next j
    Next i
    
    cnt = cnt + 1

Loop
    singleQR = a 'if you need Eigen vector ,replace a(Eigen value) to q
End Function
