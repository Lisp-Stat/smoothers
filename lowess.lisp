;;; -*- Mode: LISP; Base: 10; Syntax: ANSI-Common-Lisp; Package: SMOOTHERS -*-
(in-package #:smoothers)

;;; See https://gist.github.com/agramfort/850437
(defun lowess (x y &key
		     (f 2/3)
		     (iter 3)
		     (wt-fun (lambda (x)
			       (cube (- 1 (cube x))))))
  "Lowess smoother (robust locally weighted regression).

    Fits a nonparametric regression curve to a scatterplot.
    Parameters
    ----------
    y, x : vector
        The arrays x and y contain an equal number of elements; each pair (x[i], y[i]) defines a data point in the scatterplot.
    f : real
        The smoothing span. A larger value will result in a smoother curve.
    iter : integer
        The number of robustifying iterations. The function will run faster with a smaller number of iterations.
    Returns
    -------
    yest : vector
        The estimated (smooth) values of y."
  (let+ ((n (length x))
	 (r (ceiling (* f n)))
	 (y-est (aops:zeros n))
	 (delta (aops:ones n))
	 ((&flet clip (a min max)
	    (aops:map-array a (lambda (x) (clamp x min max)))))
	 h w)

    ;; Determine weights
    (setf h (loop
	      with vec = (make-array n :element-type 'double-float :fill-pointer 0)
	      for i from 0 below n
	      do (vector-push (aref (sort (eabs (e- x (aref x i))) #'<) r) vec)
	      finally (return vec))
	  w (aops:map-array
	     (clip (eabs (e/ (e- (aops:recycle x :outer n) (aops:recycle x :inner n))
			     (aops:recycle h :outer n)))
		   0d0 1d0)
	     wt-fun))

    (loop
      for i below iter
      with residuals and s
      do (loop
	   for k below n
	   with weights and b and A and beta
	   do (setf weights (e2* delta (select w t k))
		    b (vec 'double-float (sum (e2* weights y)) (sum (e* weights x y)))
		    A (mx 'double-float
			((sum weights) (sum (e* weights x)))
			((sum (e* weights x)) (sum (e* weights x x))))
		    beta (lla:solve A b)
			 ;; (linear-algebra:solve A b)
		    (aref y-est k) (+ (aref beta 0) (* (aref beta 1) (aref x k)))))
	 (setf residuals (e- y y-est)
	       s (median (eabs residuals))
	       delta (clip (e/ residuals (* 6 s)) -1 1)
	       delta (esquare (e- 1 (esquare delta))))
      finally (return y-est))))


