# Recurrence Plots
Software for interactively exploring Recurrence Analysis as a tool for time series analysis. Live visualization of time series and their Recurrence Plots while playing around with recurrence parameters.

The software can be used interactively and as a command line tool for time series processing. I wrote both parts for my master's thesis. Therefore, this is not an attempt to implement a comprehensive tool set for Recurrence Analysis (JRPs, automated selection of parameters, etc.) but rather an environment to try the new concepts developed in my master's thesis. For comprehensive implementation of standard recurrence analysis, very nice python packages have emerged [2,3].

# Recurrence Analysis

Recurrence Analysis is an approach to characterize multi-dimensional time series by means of self-similarity. For instance, the length, frequencies, and number of different reocurring motifs can be used to classify time series. The very basic approach is to compute the distance between all pairs of points and plot the resulting matrix. For instance, approximate repetitions in the signal then lead to diagonal lines in this Recurrence Plot. A Recurrence Plot can be characterized by a variety of features, such as average diagonal line lengths. This is then referred to as Recurrence Quantification Analysis (RQA). The RQA measures also are related to concepts in chaos theory and theoretical physics. A good introduction is given in [1].

## Deviation from traditional RQA definitions
The definition of line structures underlying the computation of Recurrence Quantification Analysis (RQA) measures differs from the traditional definition, as also explained in my blog post [5].
Line structures that are incident with the borders of the Recurrence Matrix are only considered in measures that utilize the minimum length of a line (like Determinism or Laminarity).
They are not used in measures that utilize the concrete length of a line (since the endpoint of the line lies outside the Recurrence Matrix).
However, for larger Recurrence Plots the difference should be negligible.

In addition, the software computes, visualizes, and exports a new approach to recurrence quantification,the Conditional Recurrence Plot [4].

# Implementation

The user interface relies on JavaFX and the interface code uses Java 8 features.
Time series are read in a simple CSV-based format: A d-dimensional time series of length k needs to be stored in a text file with n rows, each containing k comma-separated values.

## Maturity and Test Coverage

The computation of RQA measures has been thoroughly tested. The interfaces, however, are not at a production level, (both graphical and command line). See also the disclaimer below.

## Performance
The core RQA computations are reasonably fast, the user interface is not heavily performance optimized.
The software has some functionality to run computations in parallel and write the results into a relational database for larger scale experiments. 

# References

[1] Marwan, Norbert, M. Carmen Romano, Marco Thiel, and Jürgen Kurths. "Recurrence plots for the analysis of complex systems." Physics reports 438, no. 5-6 (2007): 237-329.

[2] https://pypi.org/project/PyRQA/

[3] http://www.pik-potsdam.de/~donges/pyunicorn/

[4] http://carl-witt.de/conditional-recurrence-time/

[5] http://carl-witt.de/edge-cases-in-rqa/

# Disclaimer

THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
