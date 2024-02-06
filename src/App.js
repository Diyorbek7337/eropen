import Container from 'react-bootstrap/Container'
import './App.css';
import Contact from './components/contact/Contact';
import Main from './components/main/Main';
import Swipers from './components/swiper/Swiper'
import {useState, useEffect} from 'react'
import Social from './components/social/Social';
import Logo from './image/logo.png'



function App() {
  const [open , setOpen] = useState(true)
  useEffect(()=> {
    setTimeout(() => {
      setOpen(false) 
    }, 2000);
  })
  const ScrollAnime = () => {
    const animItems = document.querySelectorAll(".--anim-items");
    if (animItems.length > 0) {
      function animOnScroll() {
        for (let t = 0; t < animItems.length; t++) {
          const n = animItems[t],
            i = n.offsetHeight,
            r = e(n).top,
            a = 12;
          let s = window.innerHeight - i / a;
          i > window.innerHeight && (s = window.innerHeight - window.innerHeight / a) 
          window.pageYOffset > r - s && window.pageYOffset < r + i ? n.classList.add("--animate") : n.classList.remove("--animate")
        }
  
        function e(e) {
          const t = e.getBoundingClientRect(),
            n = window.pageXOffset  || document.documentElement.scrollLeft;
          return window.scrollTo = window.pageYOffset || document.documentElement.scrollTop, {
            top: t.top + window.scrollTo,
            left: t.left + n
          }
        }
      }
      window.addEventListener("scroll", animOnScroll)
      setTimeout(() => {
        animOnScroll()
      }, 1000) 
    }
  }
  useEffect(()=>{
    ScrollAnime()
  }, [])


 
  
    return (
      <Container>
        <div className="App">
          {open ? 
            <div className="loader">
              <img src={Logo} alt="Dora Line" className='loader-image'/>
            </div>
            :
            null}
          <Main/>
          <Swipers/>
          <Social/>
          <Contact/>
        </div>
      </Container>
    );
  }

export default App
