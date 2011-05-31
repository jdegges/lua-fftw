require 'fftw'

msg = "simple interleaved in place 1d dft..."

a = {10, 0, 10, 0, 10, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}


forward_plan = fftw.plan_dft_1d (#a/2)
backward_plan = fftw.plan_dft_1d (#a/2, fftw.BACKWARD)

b = forward_plan:execute_dft (a)
c = backward_plan:execute_dft (b, 2/16)

if #a ~= #b or #b ~= #c then
  print (msg .. "failed! (1)")
  error ()
end

for i = 1, #a do
  if 1e-6 < math.abs (a[i] - c[i]) then
    print (msg .. " failed! (2)")
    error ()
  end
end

print (msg .. " success!")
