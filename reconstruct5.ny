;nyquist plug-in
;version 4
;type process
;name "Sine Peak Reconstructor (FINAL – Sample Accurate)"
;control count "Number of Sine Waves" int "" 5 1 100
;control density "Analysis Density (Snapshots)" int "" 20 1 100
;control preserve "Preserve Original Amplitude" choice "No,Yes" 0

(defun pow2-floor (n)
  (let ((p 1))
    (while (<= (* p 2) n)
      (setf p (* p 2)))
    p))

(defun hann-window (arr len)
  (dotimes (i len)
    (let ((w (* 0.5 (- 1 (cos (/ (* 2 pi i) (- len 1)))))))
      (setf (aref arr i) (* w (aref arr i))))))

(defun interpolate-peak (m1 m2 m3)
  (let ((den (- (+ m1 m3) (* 2 m2))))
    (if (= den 0)
        0.0
        (/ (- m1 m3) (* 2.0 den)))))

(defun make-zero-padded (samples winlen fftlen)
  (let ((zp (make-array fftlen)))
    (dotimes (i fftlen)
      (setf (aref zp i) 0.0))
    (dotimes (i winlen)
      (setf (aref zp i) (aref samples i)))
    zp))

(defun process-channel (sig n dens preserve)

  (let* ((sr (snd-srate sig))

         ;; EXACT sample count of selection
         (sel-samps (snd-length sig ny:all))
         (duration (/ sel-samps sr))

         ;; analysis window
         (winlen (min 65536 (pow2-floor (max 2048 sel-samps))))
         (fftlen (min 262144 (* 4 winlen)))

         (half-len (/ fftlen 2))
         (bin-width (/ sr fftlen))

         (total-mags (make-array half-len))
         (skip-dist (/ duration dens)))

    ;; initialize spectrum accumulator
    (dotimes (i half-len)
      (setf (aref total-mags i) 0.0))

    ;; ---- ANALYSIS ----
    (dotimes (s dens)
      (let* ((start-time (* s skip-dist))
             (end-time (+ start-time (/ winlen sr)))
             (chunk (extract start-time end-time (snd-copy sig)))
             (samples (snd-fetch-array chunk winlen winlen)))

        (when samples
          (hann-window samples winlen)

          (let* ((zp (make-zero-padded samples winlen fftlen))
                 (fft-data (snd-fft (snd-from-array 0.0 sr zp)
                                    fftlen fftlen nil)))

            (dotimes (i half-len)
              (let* ((re (aref fft-data (* 2 i)))
                     (im (aref fft-data (1+ (* 2 i))))
                     (mag (sqrt (+ (* re re) (* im im)))))
                (setf (aref total-mags i)
                      (+ (aref total-mags i) mag))))))))

    ;; ---- RESYNTHESIS ----
    (let ((out (make-array sel-samps)))

      ;; initialize output buffer
      (dotimes (i sel-samps)
        (setf (aref out i) 0.0))

      ;; build additive synthesis directly into array
      (dotimes (k n)

        (let ((max-val -1.0)
              (max-idx 0))

          ;; find strongest peak
          (do ((j 2 (+ j 1)))
              ((= j (- half-len 2)))
            (let ((curr (aref total-mags j)))
              (when (and (> curr max-val)
                         (> curr (aref total-mags (- j 1)))
                         (> curr (aref total-mags (+ j 1))))
                (setf max-val curr)
                (setf max-idx j))))

          (when (> max-val 0)

            (let* ((m1 (aref total-mags (- max-idx 1)))
                   (m2 (aref total-mags max-idx))
                   (m3 (aref total-mags (+ max-idx 1)))
                   (delta (interpolate-peak m1 m2 m3))
                   (true-bin (+ max-idx delta))
                   (freq (* true-bin bin-width))

                   ;; generate sine directly into buffer
                   (phase 0.0)
                   (phase-inc (* 2 pi freq (/ 1.0 sr))))

              (dotimes (i sel-samps)
                (setf (aref out i)
                      (+ (aref out i)
                         (* m2 (sin phase))))
                (setf phase (+ phase phase-inc))))

            ;; wipe neighbors
            (setf (aref total-mags max-idx) 0.0)
            (setf (aref total-mags (- max-idx 1)) 0.0)
            (setf (aref total-mags (+ max-idx 1)) 0.0))))

      ;; ---- SCALING ----
      (let* ((snd (snd-from-array 0.0 sr out))
             (out-pk (peak snd 1000000))
             (scaled
              (cond
                ((<= out-pk 0)
                 snd)

                ((= preserve 1)
                 (let ((orig-pk (peak sig 1000000)))
                   (if (> orig-pk 0)
                       (scale (/ orig-pk out-pk) snd)
                       snd)))

                (t
                 (scale (/ 0.8 out-pk) snd)))))

        ;; RETURN SAMPLE-EXACT SOUND
        scaled))))

;; ---- stereo handling ----
(if (arrayp *track*)
    (vector
      (process-channel (aref *track* 0) count density preserve)
      (process-channel (aref *track* 1) count density preserve))
    (process-channel *track* count density preserve))
