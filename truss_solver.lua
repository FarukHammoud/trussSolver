local matrix = require "matrix"
Vector = {}
function Vector:new(x,y)
	local n_vector = setmetatable({},{__index = Vector})
	n_vector.x = x
	n_vector.y = y
	return n_vector
end
function t(wik,wjk)
	if math.abs(wik) > math.abs(wjk) then
		return -wjk/wik
	else
		return -wik/wjk
	end
end
function s(wik,wjk)
	if math.abs(wik) > math.abs(wjk) then
		return c(wik,wjk)*t(wik,wjk)
	else
		return 1.0/math.sqrt(1.0+t(wik,wjk)^2)
	end
end
function c(wik,wjk)
	if math.abs(wik) > math.abs(wjk) then
		return 1.0/math.sqrt(1.0+t(wik,wjk)^2)
	else
		return s(wik,wjk)*t(wik,wjk)
	end
end
function givens(W,i,j,c,s)
	local n = #W
	local m = #W[1]
	for r = 1, m do
		local aux = c * W[i][r] - s*W[j][r]
		W[j][r] = s * W[i][r] + c * W[j][r]
		W[i][r] = aux
	end
end
function solve(W,b)
	local n = #W
	local m = #W[1]
	x = matrix:new(m,1,0)
	for k = 1, m do
		for j = n, k+1, -1 do
			local i = j - 1
			if W[j][k] ~= 0 then
				local c_ = c(W[i][k],W[j][k])
				local s_ = s(W[i][k],W[j][k])
				givens(W,i,j,c_,s_)
				givens(b,i,j,c_,s_)
			end
		end
	end
	for k = m, 1, -1 do
		local summation = 0
		for j = k+1, m do
			summation = summation + W[k][j]*x[j][1]
		end
		x[k][1] = (b[k][1] - summation)/W[k][k]
	end
	return x
end
Node = {
}
function Node:new(x,y,frame,ID)
	local n_node = setmetatable({},{__index = Node})
	n_node.x = x
	n_node.y = y
	if frame == nil then
		n_node.frame = 0
	else
		n_node.frame = frame
	end
	n_node.ID = ID
	n_node.load = Vector:new(false,false)
	n_node.displacement = Vector:new(false,false)
	return n_node
end
Bar = {
}
function Bar:new(node_r,node_s,ES)
	local n_bar = setmetatable({},{__index = Bar})
	n_bar.node_r = node_r
	n_bar.node_s = node_s
	n_bar.ES = ES
	return n_bar
end
function Bar:length()
	return math.sqrt(math.pow(self.node_r.x-self.node_s.x,2)+math.pow(self.node_r.y-self.node_s.y,2))
end
function Bar:incoming()
	local dx = self.node_s.x - self.node_r.x
	local dy = self.node_s.y - self.node_r.y
	if dx > 0 then
		return math.atan(dy/dx)
	elseif dx < 0 then
		return math.pi + math.atan(dy/dx)
	else 
		if dy >= 0 then
			return math.pi/2
		else
			return -math.pi/2
		end
	end
end
function Bar:T()
	T = matrix:new(4,4,0)
	T[1][1] = math.cos(self:incoming()-self.node_r.frame)
	T[1][2] = math.sin(self:incoming()-self.node_r.frame)
	T[2][1] = -math.sin(self:incoming()-self.node_r.frame)
	T[2][2] = math.cos(self:incoming()-self.node_r.frame)
	T[3][3] = math.cos(self:incoming()-self.node_s.frame)
	T[3][4] = math.sin(self:incoming()-self.node_s.frame)
	T[4][3] = -math.sin(self:incoming()-self.node_s.frame)
	T[4][4] = math.cos(self:incoming()-self.node_s.frame)
	return T
end
function Bar:k()
	local T_inv = matrix.invert(self:T())
	local local_kT = matrix.mul(self:local_k(),self:T())
	return matrix.mul(T_inv,local_kT)
end
function Bar:local_k()
	local K = matrix:new(4,4,0)
	K[1][1] = 1
	K[1][3] = -1
	K[3][1] = -1
	K[3][3] = 1
	K = matrix.mulnum(K,self.ES/self:length())
	return K
end
Truss = {
}
function Truss:new(path)
	local n_truss = setmetatable({},{__index = Truss})
	n_truss.nodes = {}
	n_truss.bars = {}
	n_truss:load_file(path)
	n_truss:check()
	return n_truss
end
function Truss:L(IDr,IDs)
	local L = matrix:new(4,2*self.n_nodes,0)
	L[1][2*IDr-1] = 1
	L[2][2*IDr] = 1
	L[3][2*IDs-1] = 1
	L[4][2*IDs] = 1
	return L
end
function Truss:K()
	local K = matrix:new(2*self.n_nodes,2*self.n_nodes,0)
	for i,bar in ipairs(self.bars) do
		L = self:L(bar.node_r.ID,bar.node_s.ID)
		L_t = matrix.transpose(L)
		K = matrix.add(K,matrix.mul(L_t,matrix.mul(bar:k(),L))) -- SIGNAL BAR - TRUSS HERE
	end
	return K
end
function Truss:set_U(i,value)
	if i%2 == 1 then
		self.nodes[math.ceil(i/2)].displacement.x = value
	else
		self.nodes[math.ceil(i/2)].displacement.y = value
	end
end
function Truss:set_F(i,value)
	if i%2 == 1 then
		self.nodes[math.ceil(i/2)].load.x = value
	else
		self.nodes[math.ceil(i/2)].load.y = value
	end
end
function Truss:U(i)
	if i%2 == 1 then
		return self.nodes[math.ceil(i/2)].displacement.x
	else
		return self.nodes[math.ceil(i/2)].displacement.y
	end
end
function Truss:F(i)
	if i%2 == 1 then
		return self.nodes[math.ceil(i/2)].load.x 
	else
		return self.nodes[math.ceil(i/2)].load.y 
	end
end
function Truss:order()
	order = {}

	for i,node in ipairs(self.nodes) do
		if node.load.x then
			table.insert(order,2*i-1)
		end
		if node.load.y then
			table.insert(order,2*i)
		end
	end
	cut = #order
	for i,node in ipairs(self.nodes) do
		if not node.load.x then
			table.insert(order,2*i-1)
		end
		if not node.load.y then
			table.insert(order,2*i)
		end
	end

	return order, cut
end
function Truss:U_beta()
	o, cut = self:order()
	U_beta = {}
	for i = cut+1, 2*self.n_nodes do
		table.insert(U_beta,{self:U(o[i])})
	end
	return U_beta
end
function Truss:F_alpha()
	o, cut = self:order()
	F_alpha = {}
	for i = 1, cut do
		table.insert(F_alpha,{self:F(o[i])})
	end
	return F_alpha
end
function Truss:K_sub_matrices()

	local K = self:K()
	local o, cut = self:order()
	local ordenated_K = matrix:new(2*self.n_nodes,2*self.n_nodes,0)
	for i = 1, 2*self.n_nodes do
		for j = 1, 2*self.n_nodes do
			ordenated_K[i][j] = K[o[i]][o[j]] 
		end
	end

	A = matrix.subm(ordenated_K,1,1,cut,cut)
	B = matrix.subm(ordenated_K,1,cut+1,cut,2*self.n_nodes)
	C = matrix.subm(ordenated_K,cut+1,1,2*self.n_nodes,cut)
	D = matrix.subm(ordenated_K,cut+1,cut+1,2*self.n_nodes,2*self.n_nodes)

	return A,B,C,D
end

function Truss:solve()

	A,B,C,D = self:K_sub_matrices()
	U_alpha = solve(A,matrix.mul(B,self:U_beta())-self:F_alpha())
	F_beta = matrix.mul(C,U_alpha)+matrix.mul(D,self:U_beta())
	
	o,cut = self:order()
	for i = 1, cut do
		self:set_U(o[i],-U_alpha[i][1])
	end
	for i = cut+1,2*self.n_nodes do
		self:set_F(o[i],-F_beta[i-cut][1])
	end
	print('[SOLVED] The truss is solved.')
end
function Truss:check()
	local target = 2*self.n_nodes
	for i,node in ipairs(self.nodes) do
		if not node.displacement.x then
			if not node.load.x then
				print('[UNDERDEFINED] Both disp_x and load_x arent settled for node ',i)
			else
				target = target - 1 
			end
		elseif not node.load.x then
			target = target -1
		else
			print('[OVERDEFINED] Both disp_x and load_x are settled for node ',i)
		end
		if not node.displacement.y then
			if not node.load.y then
				print('[UNDERDEFINED] Both disp_y and load_y arent settled for node ',i)
			else
				target = target - 1 
			end
		elseif not node.load.x then
			target = target -1
		else
			print('[OVERDEFINED] Both disp_y and load_y are settled for node ',i)
		end
	end
	if target == 0 then
		print('[CHECKED] The truss can be solved.')
	end
end
function Truss:load_file(path)
    local lines = {}
    for line in io.lines(path) do
    	local codes = {}
    	for item in string.gmatch(line, "%S+") do
      		table.insert(codes,item)
      	end
      	table.insert(lines,codes)
    end
    for i,line in ipairs(lines) do
    	if line[1] == 'N' then
    		ID = #self.nodes+1
    		self.nodes[ID] = Node:new(tonumber(line[2]),tonumber(line[3]),tonumber(line[4]),ID)
    	elseif line[1] == 'B' then
    		self.bars[#self.bars+1] = Bar:new(self.nodes[tonumber(line[2])],self.nodes[tonumber(line[3])],tonumber(line[4]))
    	elseif line[1] == 'F' then
    		if line[3] == 'X' then
    			self.nodes[tonumber(line[2])].load.x = tonumber(line[4])
    		elseif line[3] == 'Y' then
    			self.nodes[tonumber(line[2])].load.y = tonumber(line[4])
    		end
    	elseif line[1] == 'U' then
    		if line[3] == 'X' then
    			self.nodes[tonumber(line[2])].displacement.x = tonumber(line[4])
    		elseif line[3] == 'Y' then
    			self.nodes[tonumber(line[2])].displacement.y = tonumber(line[4])	
    		end
    	end
    end
    self.n_nodes = ID
end
--truss = Truss:new('trusses/1.trs')
--truss:solve()



