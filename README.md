This code is a work in progress and contains some known issues. Here are some of the things that will be worked on.

- Reference the article;
- Fix notation differences with the article;
- Explain which files do what, and how to use them;
- When I get a file that tries every possibility, I need to make sure to manage exceptions in case no solution exists for a specific case. I also need to check the separation distance between each thing, and discard right away solutions types (mostly the CCC ones) that require a maximum separation distance. Also, I need to manage the cases where theta is 0, 90, 180 or 270 in case an inverse trig function blows up.
- In the final code, verify in potential solutions if they overlap with C.
- Make sure to use <, >, <=, >= when appropriate.
- Maybe add a way to see if the optimal path goes through C before rejecting it or accepting it as the shortest path.
- Add a check to make sure initial conditions are outside circle C.
- Need to check sufficient condition at the start of the code, just raise a warning if it is not met.
