
/* Infosome is a javascript bioinformatics library for statistics and information theory  */

    var blosum62 = {
                '*':{'*':  1, 'A': -4, 'C': -4, 'B': -4, 'E': -4,
                      'D': -4, 'G': -4, 'F': -4, 'I': -4, 'H': -4,
                      'K': -4, 'M': -4, 'L': -4, 'N': -4, 'Q': -4,
                      'P': -4, 'S': -4, 'R': -4, 'T': -4, 'W': -4,
                      'V': -4, 'Y': -4, 'X': -4, 'Z': -4},

                 'A':{'*': -4, 'A':  4, 'C':  0, 'B': -2, 'E': -1,
                      'D': -2, 'G':  0, 'F': -2, 'I': -1, 'H': -2,
                      'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1,
                      'P': -1, 'S':  1, 'R': -1, 'T':  0, 'W': -3,
                      'V':  0, 'Y': -2, 'X':  0, 'Z': -1},

                 'C':{'*': -4, 'A':  0, 'C':  9, 'B': -3, 'E': -4,
                      'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3,
                      'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3,
                      'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2,
                      'V': -1, 'Y': -2, 'X': -2, 'Z': -3},

                 'B':{'*': -4, 'A': -2, 'C': -3, 'B':  4, 'E':  1,
                      'D':  4, 'G': -1, 'F': -3, 'I': -3, 'H':  0,
                      'K':  0, 'M': -3, 'L': -4, 'N':  3, 'Q':  0,
                      'P': -2, 'S':  0, 'R': -1, 'T': -1, 'W': -4,
                      'V': -3, 'Y': -3, 'X': -1, 'Z':  1},

                 'E':{'*': -4, 'A': -1, 'C': -4, 'B':  1, 'E':  5,
                      'D':  2, 'G': -2, 'F': -3, 'I': -3, 'H':  0,
                      'K':  1, 'M': -2, 'L': -3, 'N':  0, 'Q':  2,
                      'P': -1, 'S':  0, 'R':  0, 'T': -1, 'W': -3,
                      'V': -2, 'Y': -2, 'X': -1, 'Z':  4},

                 'D':{'*': -4, 'A': -2, 'C': -3, 'B':  4, 'E':  2,
                      'D':  6, 'G': -1, 'F': -3, 'I': -3, 'H': -1,
                      'K': -1, 'M': -3, 'L': -4, 'N':  1, 'Q':  0,
                      'P': -1, 'S':  0, 'R': -2, 'T': -1, 'W': -4,
                      'V': -3, 'Y': -3, 'X': -1, 'Z':  1},

                 'G':{'*': -4, 'A':  0, 'C': -3, 'B': -1, 'E': -2,
                      'D': -1, 'G':  6, 'F': -3, 'I': -4, 'H': -2,
                      'K': -2, 'M': -3, 'L': -4, 'N':  0, 'Q': -2,
                      'P': -2, 'S':  0, 'R': -2, 'T': -2, 'W': -2,
                      'V': -3, 'Y': -3, 'X': -1, 'Z': -2},

                 'F':{'*': -4, 'A': -2, 'C': -2, 'B': -3, 'E': -3,
                      'D': -3, 'G': -3, 'F':  6, 'I':  0, 'H': -1,
                      'K': -3, 'M':  0, 'L':  0, 'N': -3, 'Q': -3,
                      'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W':  1,
                      'V': -1, 'Y':  3, 'X': -1, 'Z': -3},

                 'I':{'*': -4, 'A': -1, 'C': -1, 'B': -3, 'E': -3,
                      'D': -3, 'G': -4, 'F':  0, 'I':  4, 'H': -3,
                      'K': -3, 'M':  1, 'L':  2, 'N': -3, 'Q': -3,
                      'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3,
                      'V':  3, 'Y': -1, 'X': -1, 'Z': -3},

                 'H':{'*': -4, 'A': -2, 'C': -3, 'B':  0, 'E':  0,
                      'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H':  8,
                      'K': -1, 'M': -2, 'L': -3, 'N':  1, 'Q':  0,
                      'P': -2, 'S': -1, 'R':  0, 'T': -2, 'W': -2,
                      'V': -3, 'Y':  2, 'X': -1, 'Z':  0},

                 'K':{'*': -4, 'A': -1, 'C': -3, 'B':  0, 'E':  1,
                      'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1,
                      'K':  5, 'M': -1, 'L': -2, 'N':  0, 'Q':  1,
                      'P': -1, 'S':  0, 'R':  2, 'T': -1, 'W': -3,
                      'V': -2, 'Y': -2, 'X': -1, 'Z':  1},

                 'M':{'*': -4, 'A': -1, 'C': -1, 'B': -3, 'E': -2,
                      'D': -3, 'G': -3, 'F':  0, 'I':  1, 'H': -2,
                      'K': -1, 'M':  5, 'L':  2, 'N': -2, 'Q':  0,
                      'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1,
                      'V':  1, 'Y': -1, 'X': -1, 'Z': -1},

                 'L':{'*': -4, 'A': -1, 'C': -1, 'B': -4, 'E': -3,
                      'D': -4, 'G': -4, 'F':  0, 'I':  2, 'H': -3,
                      'K': -2, 'M':  2, 'L':  4, 'N': -3, 'Q': -2,
                      'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2,
                      'V':  1, 'Y': -1, 'X': -1, 'Z': -3},

                 'N':{'*': -4, 'A': -2, 'C': -3, 'B':  3, 'E':  0,
                      'D':  1, 'G':  0, 'F': -3, 'I': -3, 'H':  1,
                      'K':  0, 'M': -2, 'L': -3, 'N':  6, 'Q':  0,
                      'P': -2, 'S':  1, 'R':  0, 'T':  0, 'W': -4,
                      'V': -3, 'Y': -2, 'X': -1, 'Z':  0},

                 'Q':{'*': -4, 'A': -1, 'C': -3, 'B':  0, 'E':  2,
                      'D':  0, 'G': -2, 'F': -3, 'I': -3, 'H':  0,
                      'K':  1, 'M':  0, 'L': -2, 'N':  0, 'Q':  5,
                      'P': -1, 'S':  0, 'R':  1, 'T': -1, 'W': -2,
                      'V': -2, 'Y': -1, 'X': -1, 'Z':  3},

                 'P':{'*': -4, 'A': -1, 'C': -3, 'B': -2, 'E': -1,
                      'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2,
                      'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1,
                      'P':  7, 'S': -1, 'R': -2, 'T': -1, 'W': -4,
                      'V': -2, 'Y': -3, 'X': -2, 'Z': -1},

                 'S':{'*': -4, 'A':  1, 'C': -1, 'B':  0, 'E':  0,
                      'D':  0, 'G':  0, 'F': -2, 'I': -2, 'H': -1,
                      'K':  0, 'M': -1, 'L': -2, 'N':  1, 'Q':  0,
                      'P': -1, 'S':  4, 'R': -1, 'T':  1, 'W': -3,
                      'V': -2, 'Y': -2, 'X':  0, 'Z':  0},

                 'R':{'*': -4, 'A': -1, 'C': -3, 'B': -1, 'E':  0,
                      'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H':  0,
                      'K':  2, 'M': -1, 'L': -2, 'N':  0, 'Q':  1,
                      'P': -2, 'S': -1, 'R':  5, 'T': -1, 'W': -3,
                      'V': -3, 'Y': -2, 'X': -1, 'Z':  0},

                 'T':{'*': -4, 'A':  0, 'C': -1, 'B': -1, 'E': -1,
                      'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2,
                      'K': -1, 'M': -1, 'L': -1, 'N':  0, 'Q': -1,
                      'P': -1, 'S':  1, 'R': -1, 'T':  5, 'W': -2,
                      'V':  0, 'Y': -2, 'X':  0, 'Z': -1},

                 'W':{'*': -4, 'A': -3, 'C': -2, 'B': -4, 'E': -3,
                      'D': -4, 'G': -2, 'F':  1, 'I': -3, 'H': -2,
                      'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2,
                      'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11,
                      'V': -3, 'Y':  2, 'X': -2, 'Z': -3},

                 'V':{'*': -4, 'A':  0, 'C': -1, 'B': -3, 'E': -2,
                      'D': -3, 'G': -3, 'F': -1, 'I':  3, 'H': -3,
                      'K': -2, 'M':  1, 'L':  1, 'N': -3, 'Q': -2,
                      'P': -2, 'S': -2, 'R': -3, 'T':  0, 'W': -3,
                      'V':  4, 'Y': -1, 'X': -1, 'Z': -2},

                 'Y':{'*': -4, 'A': -2, 'C': -2, 'B': -3, 'E': -2,
                      'D': -3, 'G': -3, 'F':  3, 'I': -1, 'H':  2,
                      'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1,
                      'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W':  2,
                      'V': -1, 'Y':  7, 'X': -1, 'Z': -2},

                 'X':{'*': -4, 'A':  0, 'C': -2, 'B': -1, 'E': -1,
                      'D': -1, 'G': -1, 'F': -1, 'I': -1, 'H': -1,
                      'K': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1,
                      'P': -2, 'S':  0, 'R': -1, 'T':  0, 'W': -2,
                      'V': -1, 'Y': -1, 'X': -1, 'Z': -1},

                 'Z':{'*': -4, 'A': -1, 'C': -3, 'B':  1, 'E':  4,
                      'D':  1, 'G': -2, 'F': -3, 'I': -3, 'H':  0,
                      'K':  1, 'M': -1, 'L': -3, 'N':  0, 'Q':  3,
                      'P': -1, 'S':  0, 'R':  0, 'T': -1, 'W': -3,
                      'V': -2, 'Y': -2, 'X': -1, 'Z': 4}}



  var blosum_bg_dist = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072];

    var amino_dict = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-'];


    function seq_weights(json) {

        //if(!validateAlignment(json))
         //   return;

        /*  calculate position-based sequence weights using Henikoff & Henikoff '94 method */


        var seqWeights = [], freq_count = []; // 20 aminoacids, gaps won't get in

        for (var c = 0; c <= 20; c++) {

              freq_count.push(0);

            }

        for (var c = 0; c < json.length; c++) {

              seqWeights.push(0);

            }


        for(var i=0; i < json[0].seq.length; i++) {
            

            for (var j=0; j < json.length; j++) {

                if(json[j].seq[i] != '-')
                    freq_count[amino_dict.indexOf(json[j].seq[i].toUpperCase())] += 1;

            }

            var num_observed_types = 0;

            for (var j=0; j < freq_count.length; j++) {

                if (freq_count[j] > 0)
                    num_observed_types += 1;

            }

            for (var j=0; j < json.length; j++){

                var d = freq_count[amino_dict.indexOf(json[j].seq[i].toUpperCase())] * num_observed_types;
                if (d > 0)
                    seqWeights[j] += 1 / d;

            }

            //console.log(freq_count);

        }

        //console.log(seqWeights);


        for (var w=0; w < seqWeights.length; w++) {

            seqWeights[w] /= json[0].seq.length; // need to normalize this?
        }
        //console.log(seqWeights);
        return seqWeights;
    }


    function getColumn (col_num, json) {

        col = [];

        for (var i = 0; i < json.length; i++) {

            
            col.push(json[i].seq[col_num]);

        }

        return col;

    }


    Array.prototype.sum = function () {
    for(var total = 0,l=this.length;l--;total+=this[l]); return total;
    }


   function weighted_freq_count(col, seq_weights) {

        var PSEUDOCOUNT = 0.0000001;

        /*if (seq_weights.length != col.length) {
            seq_weights = new Array(col.length);
            seq_weights.push(1);
        } */

        //console.log(seq_weights);

        var aa_num = 0;
        var freq_count = []; // don't have background distribution for gaps atm.

        for(var c = 0; c < amino_dict.length-1; c++) {
            freq_count.push(PSEUDOCOUNT);
        }

        //console.log(freq_count);


        for (var i = 0; i < amino_dict.length-1; i++ ) {

          for (var j = 0; j < col.length; j++ ) {

            if (col[j] == amino_dict[i])
              freq_count[aa_num] += 1 * seq_weights[j];  


          }

          aa_num += 1;

        }


        for (var k = 0; k < freq_count.length; k++) {

          freq_count[k] = freq_count[j] / ( seq_weights.sum() * amino_dict.length * PSEUDOCOUNT );
        }

        //console.log(freq_count);
        return freq_count;

    }


    function g_penalty(col, seq_weights) {

        /*  calculates the simple gap penalty multiplier for the column. If the 
            sequences are weighted, the gaps, when penalized, are weighted 
            accordingly.    */


            if (seq_weights.length != col.length){
              seq_weights = new Array(col.length);
              seq_weights.push(1);
            }

            var gap_sum = 0;

            for (var i = 0; i < col.length; i++) {

              if(col[i] == '-')
                gap_sup += seq_weights[i];
            }


            //console.log(gap_sum);

            return 1 - (gap_sum / seq_weights.sum());


    }


    function js_divergence(json, col, gap_penalty = 1){


      var distr = blosum_bg_dist;

      var fc = weighted_freq_count(col, seq_weights(json));


      //console.log(fc);



      if (distr.length == 20) {

        // TODO: check the background distribution for gaps and change the length of freq_count acconrdingly

      }

      // Make the distribution

      var r = [];
      for (var i = 0; i < fc.length; i++) {

        r.push(0.5 * fc[i] * 0.5 * distr[i]);

      }


      //console.log(r);

      var d = 0;

      for (var i = 0; i < fc.length; i++) {

        if (parseFloat(r[i]) != parseFloat(0.0))
          d += distr[i] * Math.log2(distr[i]/r[i]);
        else if (parseFloat(distr[i]) == parseFloat(0.0))
          d += fc[i] * Math.log2(fc[i]/r[i]);
        
        else
          d += fc[i] * Math.log2(fc[i]/r[i]) + distr[i] * Math.log2(distr[i]/r[i]);


      //console.log(d);


        d /= 2;

      }
      
      /*if (gap_penalty == 1)
        return d * g_penalty(col, seq_weights(json));
      else
        return d; */

    return d;

    }



