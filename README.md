***Japanese Tranliteration to English***

Japanese 'katakana' are special symbols used to create japanese pronunciations for foreign pronouns and technical terms. A problem of interest is to recover the original english word given the japanese pronunciation associated with a set of katakana.

One pipeline for this process is a noise channel model that assumes a given katakana instance was produced by the following process:

- An english word **E** was generated
- English word **E** was distorted into an english pronunciation **EP** via a noisy channel.
- English pronunciation **EP** was distorted to japanese pronunciation **JP** via a noisy channel.

We store the statistics for each of the processes in wfsa/wfsts.:
- The probability of generating an english word can be stored in a wfsa.
- The probability of generating an english phenome sequence given a word can be stored in a wfst.
- The probability of generating a japanese phenome sequence given an english phenome sequence can be stored in a wfst.

We can learn the statistics of english word distributions from appropriate large bodies of text. This is stored in `eword.wfsa`. We can use some naive scheme for generating english phenomes from english characters, where each transition from english character to phenome-sequence possibility is given equal probability. This is stored in `eword-epron.wfst`.

Lastly, we are interested in learning alignment statistics for english pronunciations to japanese pronunciations and storing these in `epron-jpron.wfst`. We have a set of unaligned pronunciation pairs in `epron-jpron-unsupervised.data`, and will learn the alignments statistics via an EM algorithm.

Once we have these statistics stored in WFSA/WFSTs, we use the carmel sofrware package (http://www.isi.edu/licensed-sw/carmel/) to execute the Viterbi algorithm and to output the most likely english source word that generated the japanese pronunciation.

***Running***

To learn `epron-jpron-unsupervised.wfst`, run `python run.py`.

This reads in `epron-jpron-unsupervised.data` and:
- Runs an EM algorithm to compute the probabilities for the noisy channel model (the probability of seeing a japanese phenome sequence given one english phenome)
- Produces carmel-formatted `epron-jrpon-unsupervised.wfst` from this learned model
- Calculates the probabilities of each possible alignment, given the learned channel model, for the first five phenome pairings in `epron-jpron-unsupervised.data`
- Produces `epron-jpron.alignment`
- Calculates the accuracy of produced `epron-jpron.alignment` with respect to `epron-jpron.data` (accurate alignment data)

***Translation with Carmel***
Once `epron-jpron-unsupervised.wfst` has been learned, run:

<pre>
<code>
    
</code>
</pre>

Expected Output:
<pre>
<code>
    
</code>
</pre>

***TODO***
- Add carmel commands to verify translation process.
- Clean up, modularize, document EM and parallel_data.
- Seperate out data and src files.
- Move into jupyter notebook with visualizations.
- Build hooks into carmel or write own viterbi code (slower in python...).
- Build web interface.
