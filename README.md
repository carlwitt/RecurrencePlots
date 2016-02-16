# RecurrencePlots
Compute Recurrence Plots and their RQA measures in Java. This is a relatively fast implementation I used for my master's thesis. Therefore, there are lots of possible additions (JRPs, specify Recurrence Threshold via fixed Recurrence Rate, etc.).
The user interface is not performance optimized and might perform redundant computations that lead to lower performance.

##Deviation from traditional RQA definitions
The definition of line structures underlying the computation of Recurrence Quantification Analysis (RQA) measures differs from the traditional definition.
Line structures that are incident with the borders of the Recurrence Matrix are only considered in measures that utilize the minimum length of a line (like Determinism or Laminarity).
They are not used in measures that utilize the concrete length of a line (since the endpoint of the line lies outside the Recurrence Matrix).
However, for larger Recurrence Plots the difference should be negligible.

##Maturity and Test Coverage
Most of the code is not at a production level, especially the user interface (both graphical and command line).
The computation of RQA measures has been thoroughly tested. However, a bug has been found regarding the diagonal lines that are orthogonal to the main diagonal.

##Requires Java 8
The user interface relies on JavaFX and the interface code uses Java 8 features.

##Larger Scale Experiments
I wrote some functionality to run computations in parallel and write the results into a relational database. 

##Input format
A d-dimensional time series of length k needs to be stored in a text file with n rows, each containing k comma-separated values.

##Disclaimer
THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.